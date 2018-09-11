##########################################################################################################
#                  Function about calculating the dope score  of the random forest      				 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    4/24/2011														 #
#																										 #
#																										 #
#									Revised at 4/24/2011	                         					 #
#																										 #
##########################################################################################################
#! /usr/bin/perl -w
=pod
You may freely copy and distribute this document so long as the copyright is left intact. You may freely copy and post unaltered versions of this document in HTML and Postscript formats on a web site or ftp site. Lastly, if you do something injurious or stupid
because of this document, I don't want to know about it. Unless it's amusing.
=cut
 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
  if (@ARGV != 2)
    { # @ARGV used in scalar context = number of args
          print "We use modeller 9.13, make sure to use the same version\n";
#	  print("Need modeller, so we use sysbio to run this job, and put the result in tulip server!\n");
	  print("This script tries to calculate the dope score for all models, two input parameters, one is the input folder, each subfolder is models for one target, and one output folder!\n");
      print("Use dope score, convert the range to [0,1]. \nHere we should use absolute path for each model, since calculating the score needs absolute path\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME  dir_input_target_model dir_output\n");
      print("\n********** example******\n");
          print "\n********************** ab initio for validation ************************\n";
          print "perl $0 /space2/rcrg4/CASP12_training/converted_deb_casp11 /space2/rcrg4/CASP12_training/1_calculated_scores/feature_6_dope_DB\n";
	  exit(1) ;
    }
 my $starttime = localtime();
 print "\n The time started at : $starttime.\n";
 my($input_dir)=$ARGV[0];
 my($output_dir)=$ARGV[1];
 -e $output_dir || system("mkdir $output_dir");
##########################################################################################################
#              Function about openning a directory and processing the files                				 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    12/27/2011														 #
#																										 #
##########################################################################################################
opendir(DIR, "$input_dir");
my($NUM);
my($IN,$OUT);
my($path_name,$path_matrix);
my($line);
my(@tem_split,@tem2);
my(@names,@dope_score_all,@tree_count,@for_rank);
my($index_name,$index_dope);
my($i,$j,$key,$key_rank);
my($the_name);
my($file_name);
my($read_folder,$read_folder2,$read_folder_a,$read_path);
my($write_name,$write_tree);
my(@targets);
my($target,$target_name);
my($count)=0;
my($return_val,$i);
my(@missing_folder)=();
my($missing_index)=0;
my(@score_appollo)=();
my($index_score);
my($work_dir);
my($dope_script,$pycommand);
my($dope_score);
my($min,$max);
my(%hash); # this hash table tries to score the model name and the prediction score
my @files = readdir(DIR);
foreach my $file (@files)
{
	if($file eq '.' || $file eq '..')
	{
		next;
	}
##########do something to the file###################
print "Processing $file...\n";
    $read_folder=$output_dir."/".$file.".dope_score";        # this is the output file path
    $read_folder2=$input_dir."/".$file;                      # this is the folder for each target
	if(!-s $read_folder2)
	{
		print "The target file $read_folder2 is not exists!\n";
		exit(0);
	}

	@dope_score_all=();
	@names=();
	$index_name=0;
    opendir(DIR, "$read_folder2");
    @targets = readdir(DIR);
    foreach my $target (@targets)
    {
	    if($target eq '.' || $target eq '..')
	    {
		    next;
	    }    
        $read_path=$read_folder2."/".$target;  # name for this model
print "Process $read_path...\n";
        $work_dir=$read_folder2;
        chdir($work_dir) || die "can't change to directly $work_dir\n";

system("pwd");


        $dope_script = "assess_dope.py";
        # Removes the python script file if exists
        if ( -f $dope_script)
        {
           `rm -f $work_dir/$dope_script`;
        }
        open(PY,">$dope_script") || die "can't write the script $dope_script!\n";
        print PY "from modeller import *    # Load standard Modeller classes\n";
        print PY "from modeller.scripts import complete_pdb    # Load the complete_pdb class"."\n";
        print PY "\n";
        print PY "env = environ()\n";
        print PY "env.libs.topology.read(file='\$(LIB)/top_heav.lib')\n";
        print PY "env.libs.parameters.read(file='\$(LIB)/par.lib')\n";
        print PY "\n";
        print PY "mdl = complete_pdb(env,'$read_path')\n";
        print PY "atmsel = selection(mdl.chains[0])     # Select all atoms in the first chain\n";
        print PY "atmsel = selection(mdl.chains[0])\n";
        print PY "score = atmsel.assess_dope()\n";
        close(PY);
        $pycommand = "mod9v8 ".$dope_script;
        system("$pycommand > /dev/null 2>&1");
        $dope_score = 9999;
        open(DOPE_CHECK, "$work_dir/assess_dope.log") || die "Can't open dope log file.\n";
        while(<DOPE_CHECK>)
        {
             $line = $_;
             if($line =~ /DOPE score               :/)
             {
                $dope_score = substr($line,(index($line, ":")+1),(length($line)-(index($line, ":")+1)));
                $dope_score =~ s/^\s+//; #remove leading spaces
                $dope_score =~ s/\s+$//; #remove trailing spaces
                $dope_score =~ s/ //gi;
                $dope_score =~ s/[^0-9.-]//gi;
                $dope_score =~ s/[0-9]-[0-9]//gi;
             }
        }
        close DOPE_CHECK;
        `rm $work_dir/assess_dope.log`;
        `rm $dope_script`;
        $dope_score_all[$index_name]=$dope_score;
		$names[$index_name]=$target;
		$index_name++;

	}#end of inside foreach
    if($index_name<2)
	{
		print "Only $index_name models???\n";
		exit(0);
	}
	$max=$dope_score_all[0];
	$min=$dope_score_all[0];

#################  save the score result #################
    $write_tree=$output_dir."/".$file.".dope_score";
    if(-e $write_tree)
	{
	     print "the result file | : $write_tree  ...Exists!\n"; 
	}
    else
	{ 
	     open (File, "&gt;$write_tree");
	     chmod (0777, $write_tree); 
         close (File);
    }
    $OUT = new FileHandle "> $write_tree";
    if (! defined($OUT) ) 
    {
       croak "Unable to open output file: $write_tree. Bye-bye.";
       exit(1);
    }
    for($i=0;$i<$index_name;$i++)
	{
               if($dope_score_all[$i] == 9999) {next;}
		print $OUT $names[$i]."\t".$dope_score_all[$i]."\n";
	}
	$OUT->close();
}#end foreach outside

 my $endtime = localtime();
 print  "\nThe time ended at : $endtime.\n";
