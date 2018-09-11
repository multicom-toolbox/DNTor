##########################################################################################################
#                  Function about calculating the OPUS_PSP score                          				 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    3/12/2013														 #
#																										 #
#																										 #
#									Revised at 3/12/2013	                         					 #
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
 
 sub create_list($$);
 sub process_output($$@);
  if (@ARGV != 3)
    { # @ARGV used in scalar context = number of args
	  print("This script tries to calculate OPUS_PSP score for all models, the input folder has the subfolders with the target name\n");
      print("This script only for benchmarking the normalization of OPUS_PSP score. \nHere we should use absolute path for each model, since calculating the score needs absolute path\n");
          print "Revised at 10/16/2012, the OPUS_PSP score should be 1000000 in default here, maybe this will improve the performance!\n";
     print "!!!!!!!!!!!!Attention, you have to use absolute path! And we can only use the version in sysbio server\n";

	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir_OPUS_PSP(with the file opus_psp inside this folder) dir_input_target_model dir_output\n");
      print("\n********** example******\n");

          print "perl $0 /home/rcrg4/CASP12_training/tool/OPUS/OPUS_PSP /space2/rcrg4/CASP12_training/converted_deb_casp11 /space2/rcrg4/CASP12_training/1_calculated_scores/feature_8_OPUS_DB\n";
	  exit(1) ;
    }
 my $starttime = localtime();
 print "\n The time started at : $starttime.\n";
 my($dir_exe)=$ARGV[0];

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
my($tem_path,$read_folder2,$read_folder_a,$read_path);
my($write_name,$write_tree);
my(@targets);
my($target,$target_name);
my($count)=0;
my($return_val,$i);
my(@missing_folder)=();
my($missing_index)=0;
my(@score_appollo)=();
my(@names_list);
my($index_score,$sum,$i_sum);
my($work_dir);
my($rwplus_score_ret,$rwplus_score,$list_path);
my($dope_score);
my($min,$max);
my($exe)=$dir_exe."/"."opus_psp";

my(%hash); # this hash table tries to score the model name and the energy score
my @files = readdir(DIR);
foreach my $file (@files)
{
	if($file eq '.' || $file eq '..')
	{
		next;
	}
##########do something to the file###################
print "Processing $file...\n";

	$write_tree=$output_dir."/".$file.".OPUS_PSP_score";
	
    $read_folder2=$input_pdb."/".$file;  # read the models inside
    $list_path=$output_dir."/".$file.".list";
    $tem_path=$output_dir."/".$file.".TEM";
	@names_list=();                                      # the name list
    @names_list=create_list($read_folder2,$list_path);

#print "Get name list @names_list\n";
    chdir("$dir_exe");
	$return_val=system("./opus_psp < $list_path > $tem_path");
    if($return_val !=0)
	{
		print "./opus_psp < $list_path > $tem_path fails\n";
		exit(0);
	}
	process_output($tem_path,$write_tree,@names_list);

    system("rm $tem_path");
	system("rm $list_path");

 
}#end foreach outside
 sub process_output($$@)
 {
	 my($tem_input,$output,@name_lists)=@_;
	 my($IN,$OUT,$line);
	 my(@tem_split);


#print "Get namelist: @name_lists\n";

	 $OUT=new FileHandle ">$output";
	 defined($OUT) || die "cannot open $output\n";
	 $IN=new FileHandle "$tem_input";
	 my($index)=0;
	 my($start)=0;
	 defined($IN) || die "cannot open $tem_input\n";
	 while(defined($line=<$IN>))
	 {
		 chomp($line);
		 $line=~s/\s+$//;
		 if(substr($line,0,5) eq "Input")
		 {
			 $start=1;
			 next;
		 }
		 if($start == 0)
		 {
			 next;
		 }
         @tem_split=split(/\s+/,$line);
		 if(@tem_split!=4)
		 {# this is head information
			 print "skip head : $line\n";
			 next;
		 }
		 print $OUT $name_lists[$index]."\t".$tem_split[1]."\n";
		 $index++;
	 }
     $OUT->close();
	 $IN->close();
 }
 sub create_list($$)
 {
	 my($folder,$out)=@_;
     my($OUT,$path);
	 $OUT=new FileHandle ">$out";
	 my(@name_lists)=();
	 my($index)=0;
	 defined($OUT) || die "Cannot open output $out\n";

     opendir(DIR, "$folder");
     my(@targets) = readdir(DIR);
     foreach my $target (@targets)
     {
	    if($target eq '.' || $target eq '..')
	    {
		    next;
	    }
		$path=$folder."/".$target;
		$name_lists[$index]=$target;
		$index++;
		print $OUT $path."\n";
	}
	$OUT->close();
	return @name_lists;
 }



 my $endtime = localtime();
 print  "\nThe time ended at : $endtime.\n";
