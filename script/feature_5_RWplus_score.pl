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
	  print("This script tries to calculate the RWplus score for all models, the input folder has the subfolders with the target name\n");
          print "by default, the score is 9999.\n";
	  print("You should execute the perl program like this: perl $PROGRAM_NAME RWplus_dir dir_input_target_model dir_output\n");
      print("\n********** example******\n");

print "\n********************** ab initio for validation ************************\n";
          print "perl $0  /home/tool/RWplus /space2/rcrg4/CASP12_training/converted_deb_casp11 /space2/rcrg4/CASP12_training/1_calculated_scores/feature_10_RWplus_DB\n";
  	  exit(1) ;
    }
 my $starttime = localtime();
 print "\n The time started at : $starttime.\n";
 my($tool_dir)=$ARGV[0];
 my($input_pdb)=$ARGV[1];
 my($addr_output)=$ARGV[2];
##########################################################################################################
#              Function about openning a directory and processing the files                				 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    12/27/2011														 #
#																										 #
##########################################################################################################
-s $input_pdb || die "no input folder: $input_pdb!\n";
my($NUM);
my($IN,$OUT);
my($path_name,$path_matrix);
my($line);
my(@tem_split,@tem2);
my(@names,@RW_score_all,@tree_count,@for_rank);
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
my($rwplus_score_ret,$rwplus_score);
my($dope_score);
my($min,$max);
my(%hash); # this hash table tries to score the model name and the prediction score
	$write_tree=$addr_output;

	@RW_score_all=(); # initialization
	@names=();
    $index_name=0;    # index for the total number of models
    $read_folder2=$input_pdb;  # read the models inside
    opendir(DIR, "$read_folder2");
    @targets = readdir(DIR);
    foreach my $target (@targets)
    {
	    if($target eq '.' || $target eq '..')
	    {
		    next;
	    }    
        $read_path=$read_folder2."/".$target;  # name for this model
#print "Process $read_path...\n";
########### calculate RWplus score #######################
         $work_dir=$addr_output;
         chdir "$tool_dir";
         $write_tree=$work_dir.".RWplus.out";
         open (File, "&gt;$write_tree");
         chmod (0777, $write_tree);
         close (File);
         #Execute RWplus to calculate the potential energy of each model


         $rwplus_score_ret = system("$tool_dir/calRWplus $read_path > $work_dir.RWplus.out");
         if ($rwplus_score_ret != 0)
         {
              #CleanUp();
              #die "failed to execute RWplus.\n";
         }
         $rwplus_score=9999;  #initialize
         open(RWPLUS_CHECK, "$work_dir.RWplus.out") || print "Can't open RWplus output file.\n";
         while(<RWPLUS_CHECK>)
         {
              $line = $_;
              $line =~ s/\n//;
              if($line =~ /RW potential =/)
              {
                  $rwplus_score = substr($line,(index($line, "=")+1),(index($line, "k")-(index($line, "=")+1)));
                  $rwplus_score =~ s/ //gi;
                  $rwplus_score =~ s/[^0-9.-]//gi;
                  $rwplus_score =~ s/[0-9]-[0-9]//gi;
                  
              }
         }
         close RWPLUS_CHECK;
         `rm $work_dir.RWplus.out`;
##########################################################
        $names[$index_name]=$target;
        $RW_score_all[$index_name]=$rwplus_score;
		$index_name++;

	}#end of inside foreach
    if($index_name<2)
	{
		print "Only $index_name models???\n";
		exit(0);
	}
	$max=$RW_score_all[0];
	$min=$RW_score_all[0];

#################  save the score result #################
    $write_tree=$addr_output;
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
		print $OUT $names[$i]."\t".$RW_score_all[$i]."\n";
	}
	$OUT->close();

  
 my $endtime = localtime();
 print  "\nThe time ended at : $endtime.\n";
