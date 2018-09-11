##########################################################################################################
#                  Function about calculating the RWplus score  of the random forest      				 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    9/10/2012														 #
#																										 #
#																										 #
#									Revised at 9/10/2012	                         					 #
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
  if (@ARGV != 3)
    { # @ARGV used in scalar context = number of args
	  print("This script tries to calculate the Dfire score for all models, the input folder has the subfolders with the target name, if the dfire fails, we set the score to be very bad model\n");
      print("Use dope score, convert the range to [0,1]. \nHere we should use absolute path for each model, since calculating the score needs absolute path\n");
          print "Revised at 10/16/2012, the dfire score should be 10000 in default here, maybe this will improve the performance!\n";
	  print("You should execute the perl program like this: perl $PROGRAM_NAME Dfire_dir dir_input_target_model dir_output\n");
      print("\n********** example******\n");

          print "\n********************** ab initio for validation ************************\n";
          print "perl $0 /home/tool/DFIRE2.1/ ../converted_deb_casp11 ../1_calculated_scores/feature_5_dfire_DB\n";
  	  exit(1) ;
    }
 my $starttime = localtime();
 print "\n The time started at : $starttime.\n";
 my($tool_dir)=$ARGV[0];
 my($input_pdb)=$ARGV[1];
 my($output_dir)=$ARGV[2];
##########################################################################################################
#              Function about openning a directory and processing the files                				 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    12/27/2011														 #
#																										 #
##########################################################################################################
if(!-s $output_dir)
{
	system("mkdir $output_dir");
}
-s $input_pdb || die "no input folder: $input_pdb!\n";
opendir(DIR, "$input_pdb");
my($NUM);
my($IN,$OUT);
my($path_name,$path_matrix);
my($line);
my(@tem_split,@tem2);
my(@names,@DFIRE_score_all,@tree_count,@for_rank);
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
my($index_score,$sum,$i_sum);
my($work_dir);
my($rwplus_score_ret,$rwplus_score);
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
#    $read_folder=$input_dir."/".$file."/Forest_matrix_initial";
	$write_tree=$output_dir."/".$file.".DFIRE_score";
        if(-s $write_tree)
        {# you already have this result
			print "The score is already generated, we skip here! $file \n";
            next;
        }
####################################################

	@DFIRE_score_all=(); # initialization
	@names=();
    $index_name=0;    # index for the total number of models
    $read_folder2=$input_pdb."/".$file;  # read the models inside
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
########### calculate RWplus score #######################
         $work_dir=$output_dir;
        # chdir "$tool_dir";
         $write_tree=$work_dir."/DFIRE.out";
         open (File, "&gt;$write_tree");
         chmod (0777, $write_tree);
         close (File);

         #Execute RWplus to calculate the potential energy of each model
         $rwplus_score_ret = system("$tool_dir/dfire2 $read_path > $work_dir/DFIRE.out");
         if ($rwplus_score_ret != 0)
         {

              print "failed to run $tool_dir/dfire2 $read_path > $work_dir/DFIRE.out\n";              #CleanUp();
              #die "failed to execute RWplus.\n";
         }
         $rwplus_score=9999;  #initialize
         open(RWPLUS_CHECK, "$work_dir/DFIRE.out") || print "Can't open DFIRE output file.\n";
         while(<RWPLUS_CHECK>)
         {
              $line = $_;
              $line =~ s/\n//;
			  @tem_split=split(/\s+/,$line);
			  $rwplus_score=$tem_split[1];
         }
         close RWPLUS_CHECK;
         `rm $work_dir/DFIRE.out`;
##########################################################
        $names[$index_name]=$target;
        $DFIRE_score_all[$index_name]=$rwplus_score;
		$index_name++;

	}#end of inside foreach
    if($index_name<2)
	{
		print "Only $index_name models???\n";
		exit(0);
	}
        $sum=0;
        $i_sum=0;
	$max=$DFIRE_score_all[0];
	$min=$DFIRE_score_all[0];
#################  save the score result #################
    $write_tree=$output_dir."/".$file.".DFIRE_score";
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
                if($DFIRE_score_all[$i] == 9999) {next;}          # skip it if we found the dfire fails
		print $OUT $names[$i]."\t".$DFIRE_score_all[$i]."\n";
	}
	$OUT->close();
}#end foreach outside

 my $endtime = localtime();
 print  "\nThe time ended at : $endtime.\n";
