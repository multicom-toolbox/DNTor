##########################################################################################################
#                  Function about calculating the GOAP score                            				 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    1/23/2013														 #
#																										 #
#																										 #
#									Revised at 1/23/2013	                         					 #
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
	  print("This script tries to calculate the GOAP score for all models, the input folder has the subfolders with the target name\n");
      print "Attention, the path used here must be absolute path ! Since this program's problem!\n";
	  print("You should execute the perl program like this: perl $PROGRAM_NAME addr_GOAP_tool dir_input_target_model dir_output\n");
      print("\n********** example******\n");


      print "\n******** for casp11 human prediction *********\n";
      print "perl $0 /home/tool/goap-alone/goap /space2/rcrg4/CASP12_training/converted_deb_casp11 /space2/rcrg4/CASP12_training/1_calculated_scores/feature_7_GOAP_DB\n";	  
  exit(1) ;
    }
 my $starttime = localtime();
 print "\n The time started at : $starttime.\n";
 my($addr_goap)=$ARGV[0];
 my($input_dir)=$ARGV[1];
 my($output_dir)=$ARGV[2];

 -e $input_dir || die "$input_dir not exists!\n";
 -e $addr_goap || die "$addr_goap not exists!\n";


######################################################################################################## #\
 my($IN,$OUT,$i,$line,$read_folder,$target,$write_goap,$tem_folder_root,$tem_folder,$write_path,$read_fasta,$return_val,$read_model,$path_s,$index,$copy_model);
 my(@tem_split,@targets,@goap_scores,@names,@files);


########you may change the path here #################################
 my($path_for_goap)="/home/tool/goap-alone";
######################################################################


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



#print "Input model folder is $input_dir!\n";
#print "Check succeed, start ...\n";

opendir(DIR, "$input_dir");
@files=readdir(DIR);
foreach my $file (@files)
{
	if($file eq '.' || $file eq '..')
	{
		next;
	}
    print "Processing $file...\n";
 

	$read_folder=$input_dir."/".$file;

	$write_path=$output_dir."/".$file.".TMP";                  # temporary file

    $write_goap=$output_dir."/".$file.".goap_score";         # the path for store the goap score of this target $file
	if(-s $write_goap)
	{# already processed
		print "$file already processed, since $write_goap exists! Next...\n";
		next;
	}
	else
	{
	     open (File, "&gt;$write_goap");
	     chmod (0777, $write_goap);
         close (File);		
	}


    @goap_scores=();         # put the goap scores for each model here
	@names=();               # put the target names here
	$index=0;                # index for each model


########### use goap to process each model ###############
	opendir(DIR, "$read_folder");
    @targets = readdir(DIR);
    foreach $target (@targets)
    {
	    if($target eq '.' || $target eq '..')
	    {
		    next;
	    }
		@tem_split=split(/\./,$target);
		if($tem_split[1] eq "inp")
		{# we skip the file *.inp, this is special for this program
			print "Attention, we get file $target!\n";
			#system("rm $target");
			next;
		}
        $tem_folder_root=$input_dir."/".$file."/".$target.".inp";            # temporary file for the goap score ,
	    if(-s $tem_folder_root)
	    {# already processed
		   print "$file already processed, since $tem_folder_root exists! \n";
		   #next;
	    }
	    else
	    {
	       open (File, "&gt;$tem_folder_root");
	       chmod (0777, $tem_folder_root);
           close (File);		
	    }
        $OUT = new FileHandle ">$tem_folder_root"; 
	    defined($OUT) || die "Cannot open $tem_folder_root!\n";
	    print $OUT $path_for_goap."\n";

		$read_model=$read_folder."/".$target;             # the path for the input model
		### save it to $tem_folder_root, then calculate the goap score ####
		
        print $OUT $target."\n";
        $OUT->close();

	    $goap_scores[$index]=-1;
		$names[$index]=$target;
################ use goap to predict the quality score ##################
        chdir("$read_folder");
		$tem_folder_root=$target.".inp";
        $return_val=system("$addr_goap < $tem_folder_root > $write_path");
	    if($return_val!=0)
	    {
		     print "$addr_goap < $tem_folder_root > $write_path fails, check!\n";
			 
		}
        system("rm $tem_folder_root");
        ###### read the goap result file, get the goap score for this target #######
		
        if(!-s $write_path)
		{
			$index++;
			next;
		}
        $IN=new FileHandle "$write_path";
		defined($IN) || die "$write_path cannot be opened\n";
		while(defined($line=<$IN>))
		{
			chomp($line);
			$line=~s/\s+$//;
			@tem_split=split(/\s+/,$line);
			if(@tem_split <5)
			{
				next;
			}
			$goap_scores[$index]=$tem_split[4];
		}
		$IN->close();
#########################################################################
		$index++;
	}# end of inside foreach
    
######## output the proq score for all models into this goap file ###########
    $OUT=new FileHandle">$write_goap";
	defined($OUT) || die "Cannot open output file $write_goap\n";
	for($i=0;$i<$index;$i++)
	{
		print $OUT $names[$i]."\t".$goap_scores[$i]."\n";
	}
	$OUT->close();

#	system("rm $write_path");
}#end outside for each
#system("rm -R $tem_folder");
 my $endtime = localtime();
 print  "\nThe time ended at : $endtime.\n";
