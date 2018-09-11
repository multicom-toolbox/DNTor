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
 sub read_spinex($); 

 sub read_dssp($);              # read the secondary from dssp parsed file
 sub read_spinex($);            # read the secondary structure from spineX output file
 sub cal_penalty($$);        # get the secondary structure penalty for two sequence
 sub convert_dssp_to_three_ss($); # convert a secondary structure of dssp into helix(G,H,I), strand(E,B), and loop(others).
 sub cal_normalized_penalty($$$$$$); # calculate the secondary structure penalty for two sequences with unequal sequence length, and normalize the similarity score
  if (@ARGV < 7)
    { # @ARGV used in scalar context = number of args
	  print "Check << capturing native/native like structures with a physico-chemical metric in protein folding >> \n";

	  print"This script calculates the secondary structure penalty of the secondary structure parsed by DSSP and secondary structure predicted by spineX for one target, if you need to process more than one target, you need to call this script several times. Use the percentage of secondary structure similarity as the quality score!\n";
	  print "\n************** Renzhi Cao *******************\n";
	  print "Input:\n";
	  print "0. Dir of all models (for one target)!\n";
	  print "1. address of fasta sequence \n";
	  print "2. address of LCS \n";
	  print "3. address of dsspcmbi\n";
	  print "4. address of dssp2dataset.pl \n";
	  print "5. address of spX.pl\n";
	  print "6. Dir of output\n";
          print "7. Optional. one, only for native pdb with one folder\n";      

	  print "\n**************** Renzhi Cao *****************\n";
	  print "\nFor example:\n";
	  print "perl $PROGRAM_NAME ../data/casp10_stage1_server_prediction_pdb/T0644 /exports/store1/rcrg4/REAL_NATIVE_from_CASP/data_downloaded_from_CASP/data_sequence/casp10_seq/T0644.fasta /exports/store1/rcrg4/tool_I_develop/LCS /exports/store1/tool/dsspcmbi ../Feature_generate_script/dssp2dataset.pl /exports/store1/tool/spine_X/spX.pl ../test/feature_2_2_test\n";
	  exit(0);
	}
############### the copy right for each script #############################
print  "REM       \n";
print  "REM          *********************************************************************************\n";
print  "REM          *                                                                               *\n";
print  "REM          *                                                                               *\n";
print  "REM          *      Developed By :   Renzhi Cao (rcrg4\@mail.missouri.edu)                   *\n";
print  "REM          *      Copyright    :   Dr. Jianlin Cheng's BDML Laboratory                     *\n";
print  "REM          *      Release Date :   Feb 2, 2014                                             *\n";
print  "REM          *      Vesion       :   1.0                                                     *\n";
print  "REM          *                                                                               *\n";
print  "REM          *********************************************************************************\n";
print  "REM       \n";
############################################################################

 my($dir_models)=$ARGV[0];
 my($addr_seq)=$ARGV[1];
 my($addr_LCS)=$ARGV[2];
 my($addr_dssp)=$ARGV[3];
 my($addr_dssp_perl)=$ARGV[4];
 my($addr_spx)=$ARGV[5];
 my($dir_output)=$ARGV[6];

 -s $dir_models || die "cannot open input $dir_models\n";
 -s $addr_seq || die "cannot open input $addr_seq\n";
 -s $addr_LCS || die "cannot open input $addr_LCS\n";
 -s $addr_dssp || die "cannot open input $addr_dssp\n";
 -s $addr_dssp_perl || die "cannot open input $addr_dssp_perl\n";
 -s $addr_spx || die "cannot open input $addr_spx\n"; 
 -s $dir_output || system("mkdir $dir_output");


 my($one_model) = 0;
 my($tmp_targets) = $dir_output.".tmp";
 system("mkdir $tmp_targets");
 if(@ARGV>7)
 {
      if($ARGV[7] eq "one")
      {
          print "great, only acess one model\n";
          $one_model= 1;
          system("cp $dir_models $tmp_targets/");
          $dir_models = $tmp_targets;
      }
 }



 my($path_target,$sequence_name,$spx_seq_dir,$return_val);
 my(@tem_split,@targets,@files,@searches);
 my($IN,$line,$i,$OUT);

 @tem_split=split(/\//,$addr_seq);
 $sequence_name=$tem_split[@tem_split-1];         # get the sequence name

 my($target_name);
 @tem_split=split(/\//,$dir_models);
 $target_name=$tem_split[@tem_split-1];         # get the target name
######## use spineX to predict the secondary structure for the input sequence ############
 # check already finished 
 my($check_path) = $dir_output."/"."spineX_out"."/".$sequence_name.".spXout";
 if(-s $check_path)
 {
    print "already processed the output, skip!\n";
    exit(0);
#    return 1;
 } 
 
 my($tem_output)=$dir_output."/"."spineX_out";
 
 
 $spx_seq_dir=$dir_output."/"."sequence_folder";
 -s $spx_seq_dir || system("mkdir $spx_seq_dir");

 $return_val = system("cp $addr_seq $spx_seq_dir/");
 if($return_val !=0)
 {
         print "cp $addr_seq $spx_seq_dir/ fails\n"; 
         exit(0);
 }

 my($tem_list)=$dir_output."/"."spineX_list";
 $OUT = new FileHandle ">$tem_list";
 print $OUT $sequence_name."\n";
 $OUT->close();
 


 $return_val = system("perl $addr_spx $tem_list $spx_seq_dir $tem_output");
 if($return_val !=0)
 {
         print "perl $addr_spx $tem_list $spx_seq_dir $tem_output fails\n";

 }


 my($spinex_ss_file)=$tem_output."/".$sequence_name.".spXout";
 if(!-s $spinex_ss_file)
 {
	 print "$spinex_ss_file not exists, check why spineX fails! \n";
	 exit(0);
 }
 my($spi_seq,$spi_ss)=read_spinex($spinex_ss_file);
 
# print "You get \n$spi_seq\n$spi_ss\n";

 ####### 2. Now we use dssp to process each model, and get the parsed secondary structure for each model #########
 my($path_write,$write_dssp,$dssp_parsed,$target);
 my($tem_dssp_out)=$dir_output."/"."dssp_output";
 -s $tem_dssp_out || system("mkdir $tem_dssp_out");
 
opendir(DIR, "$dir_models");
@targets=readdir(DIR);
foreach $target (@targets)
{
	if($target eq '.' || $target eq '..')
	{
		next;
	}
    print "Processing $target...\n";
	$path_target=$dir_models."/".$target;
    $write_dssp=$tem_dssp_out."/".$target.".dsspout";         # the path for store the dssp output of this target $file
	if(-s $write_dssp)
	{# already processed
		#print "$target already processed, since $write_dssp exists! next ...\n";
		#next;
	}
	else
	{
	     open (File, "&gt;$write_dssp");
		 chmod (0777, $write_dssp);
         close (File);		
	}
    $return_val=system("$addr_dssp $path_target $write_dssp");
	if($return_val!=0)
	{
		print "$addr_dssp $path_target $write_dssp fails!\n";
	    #exit(0);
		next;
	}	     
	$dssp_parsed=$tem_dssp_out."/".$target.".dssp_parsed";  # the path for parsed dssp output result
	if(-s $dssp_parsed)
	{# already processed
		#print "$file already processed, since $dssp_parsed exists! next ...\n";
		#next;
	}
	else
	{
		 open (File, "&gt;$dssp_parsed");
	     chmod (0777, $dssp_parsed);
	     close (File);		
	}
	$return_val=system("perl $addr_dssp_perl $write_dssp $dssp_parsed");
	if($return_val!=0)
	{
		print "perl $addr_dssp_perl $write_dssp $dssp_parsed fails!\n";
		#exit(0);
		next;
	}
	system("rm $write_dssp");           # remove the original dssp output, we can just use the parsed dssp output .
}# end of outside foreach
   
  
 ######## calculate the secondary structure similarity between spineX predicted and dssp parsed #########

my($dssp_aa,$dssp_ss,$similarity_score,$path_LCS,$path_real,$path_parse);
my($LCS_dir)=$dir_output."/"."LCS";
-s $LCS_dir || system("mkdir $LCS_dir");

my($path_ss_similarity)=$dir_output."/".$target_name.".ss_similarity";
my($OUT_SS);
$OUT = new FileHandle ">$path_ss_similarity";
opendir(DIR, "$dir_models");
@targets=readdir(DIR);
foreach $target (@targets)
{
	if($target eq '.' || $target eq '..')
	{
		next;
	}
    $path_target = $tem_dssp_out."/".$target.".dssp_parsed";         # the dssp parsed output
	if(!-s $path_target)
	{
		print "not existing $path_target!\n";
		next;
	}
    ($dssp_aa,$dssp_ss)=read_dssp($path_target);      # get the amino acids and secondary struture from the dssp parsed output
    
    
	#print "We get dssp parsed for $target:\n";
	#print "$dssp_aa\n$dssp_ss\n";
    if($spi_seq eq $dssp_aa)
	{
		$similarity_score = cal_penalty($dssp_ss,$spi_ss);        # calculate the secondary penalty score, no normalization!!!
	}
	else
	{# we need to align the sequence and then normalize the similarity score based on the spineX predicted secondary structure
        $path_LCS=$LCS_dir."/".$target.".alignment";
        $path_real = $LCS_dir."/".$target.".real";
		$path_parse = $LCS_dir."/".$target.".parsed";
		
        $OUT_SS = new FileHandle ">$path_real";
		print $OUT_SS ">$target released sequence\n";
		print $OUT_SS $spi_seq."\n";
		$OUT_SS->close();
        $OUT_SS = new FileHandle ">$path_parse";
		print $OUT_SS ">$target parsed from model sequence\n";
		print $OUT_SS $dssp_aa."\n";
		$OUT_SS->close();		
		$similarity_score = cal_normalized_penalty($addr_LCS,$path_real,$path_parse,$path_LCS,$spi_ss,$dssp_ss);
	}
    
	if($similarity_score == -1)
	{# the program fails, return 
		next;
	}
    print $OUT $target."\t".$similarity_score."\n";

}
$OUT->close();
system("rm -R $tmp_targets");


 sub cal_normalized_penalty($$$$$$)
 {
	 my($addr_LCS,$path_real,$path_parse,$path_LCS,$real_ss,$dssp_ss)=@_;
	 my($return_val)=system("$addr_LCS $path_real $path_parse > $path_LCS");
	 if(!-s $path_LCS)
	 {
		 print "$addr_LCS $path_real $path_parse > $path_LCS fails!\n";
		 return -1;
	 }
	 my($total_length)=length($real_ss);      # the total length of relased sequence 
	 if($total_length == 0)
	 {
		 print "Warning, check secondary structure $real_ss, 0?\n";
		 return -1;
	 }
     my($aligned_length);
     my($IN,$line);
	 my(@tem_split);
	 my(@sequences)=();
	 my($index)=0;
	 $IN = new FileHandle "$path_LCS";
	 while(defined($line=<$IN>))
	 {
		 chomp($line);
		 if(substr($line,0,1) eq "#")
		 {
			 next;
		 }
		 $sequences[$index++]=$line;
	 }
	 $IN->close();
     if($index<4)
	 {
		 print "The output of LCS is abnormal, check $addr_LCS $path_real $path_parse > $path_LCS\n";
		 return -1;
	 }
	 
	 my($aligned_real_ss)="NULL";
	 my($aligned_dssp_ss)="NULL";
	 my($i);
	 for($i=0;$i<length($real_ss);$i++)
	 {
		 if(substr($sequences[1],$i,1) eq "-")
		 {
			 next;
		 }
		 if($aligned_real_ss eq "NULL")
		 {
			 $aligned_real_ss = substr($real_ss,$i,1);
		 }
		 else
		 {
             $aligned_real_ss.=substr($real_ss,$i,1);
		 }
	 }
	 for($i=0;$i<length($dssp_ss);$i++)
	 {
		 if(substr($sequences[3],$i,1) eq "-")
		 {
			 next;
		 }
		 if($aligned_dssp_ss eq "NULL")
		 {
			 $aligned_dssp_ss = substr($dssp_ss,$i,1);
		 }
		 else
		 {
             $aligned_dssp_ss.=substr($dssp_ss,$i,1);
		 }
	 }


	 $aligned_length = length($aligned_dssp_ss);



     $i=cal_penalty($aligned_real_ss,$aligned_dssp_ss);        # get the penalty score for the aligned part
     my($simi_score)= $i * ($aligned_length/$total_length);
	 
     return $simi_score;
 }




 sub read_spinex($)
 {
         my($input)=@_;
         my($IN,$line);
         my(@tem_split);
         my($seq,$ss);

         $IN=new FileHandle "$input";
         if(defined($line=<$IN>))
         {# skip the head information
         }
         if(defined($line=<$IN>))
         {
                 chomp($line);
                 $line=~s/\s+$//;
                 @tem_split=split(/\s+/,$line);
                 if(@tem_split<3)
                 {
                         next;
                 }
                 $ss=$tem_split[3];          # get the first ss
				 $seq=$tem_split[2];          # get the first amino acids
         }
         while(defined($line=<$IN>))
         {
                 chomp($line);
                 $line=~s/\s+$//;
                 @tem_split=split(/\s+/,$line);
                 if(@tem_split<3)
                 {
                         next;
                 }
                 $ss.=$tem_split[3];          # get the ss
				 $seq.=$tem_split[2];         # get the amino acid

     }
         $IN->close();
         return ($seq,$ss);
 }






 sub convert_dssp_to_three_ss($)
 {
	 my($cha)=@_;
	 if($cha eq "G" || $cha eq "H" || $cha eq "I")
	 {# helix
		 return "H";
	 }
	 elsif($cha eq "E" || $cha eq "B")
	 {#strand
		 return "E";
	 }
	 else
	 {#coil
		 return "C";
	 }
 }

 sub read_dssp($)
 {
	 my($input)=@_;
	 my($IN,$line);
	 $IN=new FileHandle "$input";
	 my(@aa);
	 my($i);
	 my($dssp_aa)="NULL";
	 if(defined($line=<$IN>))
	 {
		 # this is for chain ID
	 }
	 if(defined($line=<$IN>))
	 {
		 # this is for total number
	 }
	 if(defined($line=<$IN>))
	 {
		 chomp($line);
         @aa=split(/\s+/,$line);
		 if(@aa<1)
		 {# empty amino acid??? 
			 next;
		 }
		 $dssp_aa=$aa[0];
		 for($i=1;$i<@aa;$i++)
		 {
			 $dssp_aa.=$aa[$i];
		 }
		 # this is for amino acid
	 }
	 if(defined($line=<$IN>))
	 {
		 # this is for index
	 }
	 if(defined($line=<$IN>))
	 {
		 # this is secondary structure
		 chomp($line);
	 }

	 $IN->close();
	 my(@ss)=split(/\s+/,$line);
	 my($dssp_ss)="NULL";
	 if(@ss < 1)
	 {
		 return ($dssp_aa,$dssp_ss);
	 }
	 
	 $dssp_ss=convert_dssp_to_three_ss($ss[0]);
	 for($i=1;$i<@ss;$i++)
	 {
		 $dssp_ss.=convert_dssp_to_three_ss($ss[$i]);
	 }
	 return ($dssp_aa,$dssp_ss);

 }
 sub cal_penalty($$)
 {
	 my($seq1,$seq2)=@_;
	 my($i,$len,$penalty);
	 $penalty=0;

	 $len=length($seq1);
	 ###### 1. Check penalty for helix #######
	 for($i=0;$i<$len;$i++)
	 {
		 if((substr($seq1,$i,1) eq "H") && (substr($seq2,$i,1) ne "H"))
		 {
			 $penalty++;
		 }
		 elsif((substr($seq1,$i,1) ne "H") && (substr($seq2,$i,1) eq "H"))
		 {
			 $penalty++;
		 }
	 }
	 ###### 2. Check penalty for sheet #######
	 for($i=0;$i<$len;$i++)
	 {
		 if((substr($seq1,$i,1) eq "E") && (substr($seq2,$i,1) ne "E"))
		 {
			 $penalty++;
		 }
		 elsif((substr($seq1,$i,1) ne "E") && (substr($seq2,$i,1) eq "E"))
		 {
			 $penalty++;
		 }
	 }
	 ######## convert it into a score ########
	 $penalty=(2*$len-$penalty)/(2*$len);

	 return $penalty;
 }
 
=pod
 

 my(%hash)=();                                            
 my($path_target,$search_name,$path_model1,$path_model2,$path_search,$search,$path_out,$name,$path_write,$file,$return_val,$seq1,$seq2,$quality);
 my(@tem_split,@targets,@files,@searches);
 my($IN,$line,$i,$OUT);
###########################################################################################################

 my($target);
 ############# read the dssp parsed secondary structure for each model ###############
 opendir(DIR,"$addr_input1");
 @files=readdir(DIR);
 foreach $file (@files)
 {

	if($file eq '.' || $file eq '..')
	{
	  next;
    }
	$path_target=$addr_input1."/".$file;          # the target folder
print "Processing $path_target ...\n";
	
    $path_out = $dir_output."/".$file.".secondary_structure_difference";      # the output file
	$OUT=new FileHandle ">$path_out";
	defined($OUT) || die "Cannot open output $path_out\n";
    opendir(DIR,"$path_target");
	@targets=readdir(DIR);
	foreach $target (@targets)
	{
		 if($target eq '.' || $target eq '..')
		 {
			 next;
		 }
		 $path_model1=$path_target."/".$target;             # this is the dssp parsed secondary structure. 
         @tem_split=split(/\./,$target);
		 $name=$tem_split[0];                             # this is the model name
		 for($i=1;$i<@tem_split-1;$i++)
		 {
			 $name.=".".$tem_split[$i];
		 }
		 $path_search=$addr_input2."/".$file;              # search in spineX secondary structure 
		 $path_model2 = "NULL";
		 opendir(DIR,"$path_search");
		 @searches=readdir(DIR);
		 foreach $search (@searches)
		 {
			 if($search eq '.' || $search eq '..')
			 {
				 next;
			 }
			 @tem_split=split(/\./,$search);
			 $search_name=$tem_split[0];
			 for($i=1;$i<@tem_split-1;$i++)
			 {
				 $search_name.=".".$tem_split[$i];
			 }
			 if($search_name eq $name)
			 {
				 $path_model2=$path_search."/".$search;       # find the spineX secondary structure for the same model
			 }
		 }
		 if($path_model2 eq "NULL")
		 {
			 print "Why the model $path_model1 don't find the other model $path_search/$name ???\n";
			 #exit(0);
			 next;
		 }
		 ########### now read the two secondary structure and compare them ###########
		 $seq1=read_dssp($path_model1);      # read dssp result
		 
		 $seq2=read_spinex($path_model2);    # read spinex result
		 if($seq1 eq "NULL" || $seq2 eq "NULL")
		 {
			 print "The secondary structure for dssp processed $path_model1 is $seq1, spinex processed $path_model2 is $seq2, why empty??? \n";
			 next;
			 #exit(0);
		 }
         
		 #print "Sequence 1 : \n$seq1\n Sequence 2 : \n$seq2\n";
		 ######### here can be changed! For quality calculation     ###########
		 $quality = cal_difference($seq1,$seq2);           # get the difference as the quality
		 print $OUT $name."\t".$quality."\n";
	}
	$OUT->close();
 }



 sub read_spinex($)
 {
	 my($input)=@_;
	 my($IN,$line);
	 my(@tem_split);
	 my($seq);

	 $IN=new FileHandle "$input";
	 if(defined($line=<$IN>))
	 {# skip the head information
	 }
	 if(defined($line=<$IN>))
	 {
		 chomp($line);
		 $line=~s/\s+$//;
		 @tem_split=split(/\s+/,$line);
		 if(@tem_split<3)
		 {
			 next;
		 }
		 $seq=$tem_split[3];          # get the first ss
	 }
	 while(defined($line=<$IN>))
	 {
		 chomp($line);
		 $line=~s/\s+$//;
		 @tem_split=split(/\s+/,$line);
		 if(@tem_split<3)
		 {
			 next;
		 }
		 $seq.=$tem_split[3];          # get the ss

     }
	 $IN->close();
	 return $seq;
 }



