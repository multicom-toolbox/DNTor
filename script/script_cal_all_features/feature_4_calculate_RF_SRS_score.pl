##########################################################################################################
#                  Function about calculating the RF_SRS score                          				 #
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
  if (@ARGV != 4)
    { # @ARGV used in scalar context = number of args
	  print("This script tries to calculate RF_SRS score for all models, the input folder has the subfolders with the target name\n");
      print("This script only for benchmarking the normalization of RF_SRS score. \nHere we should use absolute path for each model, since calculating the score needs absolute path\n");
          print "Revised at 10/16/2012, the RF_SRS score should be 1000000 in default here, maybe this will improve the performance!\n";
	  print "\n*********************\nThe RF_SRS has two different potential file, RF_HA_SRS and RF_CB_SRS_OD\n***************\n";

	  print("You should execute the perl program like this: perl $PROGRAM_NAME exe_calc_energy addr_potential_file dir_input_target_model dir_output\n");
      print("\n********** example******\n");
  
          print "\n********************** ab initio for validation ************************\n";   
          print "perl $0 /space2/rcrg4/transfer/tool/Energy_function_RF_SRS/calc_energy /space2/rcrg4/transfer/tool/Energy_function_RF_SRS/RF_CB_SRS_OD ../converted_deb_casp11 ../1_calculated_scores/feature_4_RF_SRS_DB\n";
  	  exit(1) ;
    }
 my $starttime = localtime();
 print "\n The time started at : $starttime.\n";
 my($exe)=$ARGV[0];
 my($addr_potential)=$ARGV[1];
 my($input_pdb)=$ARGV[2];
 my($output_dir)=$ARGV[3];
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
	$write_tree=$output_dir."/".$file.".energy_score";
        if(-s $write_tree)
        {# you already have this result
			print "The score is already generated, we skip here! $file \n";
            next;
        }
=pod
	$read_path=$read_folder."/".$file.".name";
	if(!-s $read_path)
	{
		print "not exists $read_path !\n";
		exit(0);
	}
###############read model names ####################
    $IN = new FileHandle "$read_path";
    if (! defined($IN)) 
    {
       print("Can't open spec file $read_path: $!\n");
       return 0;
    }
	@names=();
	$index_name=0;
    while ( defined($line = <$IN>))
    {#get the name list
	  #read something here
	  chomp($line);
	  $line=~s/\s+$//;  # remove the windows character
	  @tem_split=split(/\s+/,$line);
	  if(@tem_split<1)
	  {
		  next;
	  }
	  $names[$index_name]=$line;
	  $index_name++;
	}
=cut
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
         
         $write_tree=$work_dir."/energy.out";
         open (File, "&gt;$write_tree");
         chmod (0777, $write_tree);
         close (File);
         #Execute RWplus to calculate the potential energy of each model
         $rwplus_score_ret = system("$exe $read_path $addr_potential > $work_dir/energy.out");
         if ($rwplus_score_ret != 0)
         {
              print "Fail to excute ./$exe $read_path $addr_potential > $work_dir/energy.out \n";
              #CleanUp();
              #die "failed to execute RWplus.\n";
         }
         $rwplus_score=9999;  #initialize
         open(RWPLUS_CHECK, "$work_dir/energy.out") || print "Can't open energy output file $work_dir/energy.out.\n";
         while(<RWPLUS_CHECK>)
         {
              $line = $_;
              $line =~ s/\n//;
			  @tem_split=split(/\s+/,$line);

			  if(@tem_split>1)
			  {
				  next;
			  }
			  $rwplus_score=$line;
         }
         close RWPLUS_CHECK;
         `rm $work_dir/energy.out`;
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
=pod
        $sum=0;
        $i_sum=0;
	$max=$DFIRE_score_all[0];
	$min=$DFIRE_score_all[0];
	for($i=0;$i<$index_name;$i++)
	{
                if($DFIRE_score_all[$i] == 10000)
                {# if the score is the initialize score
                   next; 
                }
		if($DFIRE_score_all[$i]>$max)
		{
			$max=$DFIRE_score_all[$i];
		}
		if($DFIRE_score_all[$i]<$min)
		{
			$min=$DFIRE_score_all[$i];
		}
                $sum+=$DFIRE_score_all[$i];
                $i_sum++; 
	}
        if($i_sum!=0)
        {
               $sum/=$i_sum;
        }
        for($i=0;$i<$index_name;$i++)
        {
           if($DFIRE_score_all[$i] == 10000)
           {
              $DFIRE_score_all[$i]=$max;  # use the max score as the default
           }
        }
	if($max!=$min)
	{
	    for($i=0;$i<$index_name;$i++)
	    {
	    	$DFIRE_score_all[$i]-=$min;
	    	$DFIRE_score_all[$i]/=$max-$min;
            $DFIRE_score_all[$i]=1-$DFIRE_score_all[$i];
	    }
	}
	else
	{
		$DFIRE_score_all[$i]=0.5;
	}
=cut
#################  save the score result #################
    $write_tree=$output_dir."/".$file.".energy_score";
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
                if($DFIRE_score_all[$i] == 9999) {next;}
		print $OUT $names[$i]."\t".$DFIRE_score_all[$i]."\n";
	}
	$OUT->close();
}#end foreach outside

 my $endtime = localtime();
 print  "\nThe time ended at : $endtime.\n";
