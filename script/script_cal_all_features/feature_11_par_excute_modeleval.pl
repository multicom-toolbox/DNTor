##########################################################################################################
#                             Function about using tool modeleval to get the score for all models        #
#input: input_dir(with all models in a directory for each target), input_seq_dir addr_beta_contact_map.sh#
#                                addr_model_eval.sh output_dir		                                     #
#										Renzhi Cao  													 #
#																										 #
#									    1/24/2012														 #
#																										 #
#																										 #
#									Revised at 1/25/2012												 #
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
  if (@ARGV < 5)
    { # @ARGV used in scalar context = number of args
	  print("This program using tool modeleval to get the score for all models. Can be run in parallel!\n");
      print("You should execute the perl program like this: perl $PROGRAM_NAME  input_dir(with all models in a directory for each target) input_seq_dir addr_beta_contact_map.sh addr_model_eval.sh output_dir native?(option)!\n");
      print "\n********** casp11 human prediction test **************\n";
#      print "perl $0 /exports/store1/rcrg4/CASP_new_score_prediction/data/casp10_human_all_targets /exports/store1/rcrg4/REAL_NATIVE_from_CASP/data_downloaded_from_CASP/data_sequence/casp10_seq /exports/store1/tool/betacon/bin/beta_contact_map.sh /exports/store1/tool/model_eva1.0/bin/model_eval.sh /exports/store1/rcrg4/CASP_new_score_prediction/result_for_human_prediction/casp10_human_all_modeleva_score\n";
print "\n********************** ab initio for validation ************************\n";
      print "perl $0 /space2/rcrg4/CASP12_training/converted_deb_casp11 ../sequences /home/tool/betacon/bin/beta_contact_map.sh /home/tool/model_eva1.0/bin/model_eval.sh ../1_calculated_scores/feature_11_modeleva_DB\n";
      exit(1) ;
    }
 my $starttime = localtime();
 print "\nThe time started at : $starttime.\n";
 my($input_target)=$ARGV[0];
 my($input_seq)=$ARGV[1];
 my($addr_beta)=$ARGV[2];
 my($addr_model)=$ARGV[3];
 my($output)=$ARGV[4];
 if(!-s $input_target || !-s $output || !-s $input_seq)
 {
	print "The input file $input_target or $output or $input_seq does not exists, please check!\n";
	exit(1);
 } 
#############################################################
 my($file,$target);
 my(@files,@ttt,@targets);
 my($IN_LIST);
 my($seq_path,$write_folder,$path_read,$target_path);
 my($file_name);
 my($is_native) = 0;
 if(@ARGV>5)
 {
    if($ARGV[5] eq "native")
    {
       $is_native = 1;
    }
    else
    {
       die " You must give the native as option, only for trainning data\n";
    }
 }
#############################################################
 opendir(DIR, "$input_target");
 @files = readdir(DIR);
 foreach $file (@files)
 {
	if($file eq '.' || $file eq '..')
	{
		next;
	}
	$seq_path=$input_seq."/".$file.".fasta"; # the fasta sequence input
        @ttt = split(/\./,$file);
        
	if(!-s $seq_path)
	{
                $seq_path = $input_seq."/".$ttt[0].".fasta";
                if(!-s $seq_path) 
                { 
                    print "$seq_path not exists, which contains the sequence of the target, please check!\n";   
		    next;
                }
                
	}
    $write_folder=$output."/".$file;
        if(-s $write_folder)
        {#already have the result!
            next;
        }
	`mkdir $write_folder`;   # create write folder
    $path_read=$input_target."/".$file; # the read folder
    `$addr_beta $seq_path $write_folder`;
     if($is_native == 1)
     {
          $file_name = $output."/".$file.".model_eval_score";
          `$addr_model $seq_path $write_folder $path_read > $file_name`;    # get the score for each model
          next;
     }

    opendir(DIR, "$path_read");
    @targets = readdir(DIR);
	foreach $target (@targets) 
	{
	   if($target eq '.' || $target eq '..')
	   {
     		next;
       }
######### create file to save the score #################
     $file_name=$write_folder."/".$target.".model_eval_score";
     if(-e $file_name)
		  {
		     print "the result file | : $file_name  ...Exists!\n"; 
	      }
     else 
		  { 
			 open (File, "&gt;$file_name");
		     chmod (0777, $file_name); 
             close (File);
          }
       $target_path=$path_read."/".$target;  # the target path
	
	   `$addr_model $seq_path $write_folder $target_path > $file_name`; # get the score for each model
	}
 }
 my $endtime = localtime();
 print  "\n$PROGRAM_NAME ended at : $endtime.\n";
