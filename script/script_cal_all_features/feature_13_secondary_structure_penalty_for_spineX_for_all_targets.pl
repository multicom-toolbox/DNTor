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
 sub check_it($)
 {
     my($file) = shift;
     if(!-s $file)
     {
        return 0;
     }
     my($IN) = new FileHandle "$file";
     my($line);
     my(@tem);
     if(defined($line=<$IN>))
     {
          @tem = split(/\s+/,$line);
          if(@tem<1) {return 0;}
     }
     $IN->close();
     return 1;

 }

  if (@ARGV < 8)
    { # @ARGV used in scalar context = number of args
	  
	  print"This script calculate the secondary structure penalty score parsed by DSSP and secondary structure predicted by spineX for all targets. Use the  secondary structure penalty score of helix and sheet as the quality score!\n";
	  print "\n************** Renzhi Cao *******************\n";
	  print "Input:\n";

	  print "0. Dir of all targets (each subfolder is a target tarball)!\n";
	  print "1. dir of fasta sequences \n";
	  print "2. address of feature_2_2_secondary_structure_penalty_for_spineX_for_one_target.pl \n";
	  print "3. address of LCS \n";
	  print "4. address of dsspcmbi\n";
	  print "5. address of dssp2dataset.pl \n";
	  print "6. address of spX.pl\n";
	  print "7. Dir of output\n";
           print "Optional, 8. one. for native pdbs\n";

	  print "\n**************** Renzhi Cao *****************\n";
	  print "\nFor example:\n";
print "\n********************** CASP11 for validation ************************\n";
          print "perl $0 ../converted_deb_casp11 ../sequences feature_13_secondary_structure_penalty_for_spineX_for_one_target.pl /space2/rcrg4/transfer/tool/LCS /space2/rcrg4/transfer/tool/dsspcmbi dssp2dataset.pl /home/tool/spine_X/spX.pl ../1_calculated_scores/feature_13_sp_DB\n";


	  exit(0);
	}

 
 my($dir_targets)=$ARGV[0];
 my($dir_fasta)=$ARGV[1];
 my($addr_feature_pl)=$ARGV[2];
 my($addr_LCS)=$ARGV[3];
 my($dssp)=$ARGV[4];
 my($dssp_pl)=$ARGV[5];
 my($spx)=$ARGV[6];
 my($dir_output)=$ARGV[7];

 -s $dir_output || system("mkdir $dir_output");

 my($tem_output)=$dir_output."/TEM";
 -s $tem_output || system("mkdir $tem_output");

 my($one_model) = 0;
 if(@ARGV>8)
 {
     if($ARGV[8] eq "one")
     {
         print "setting one model ...\n";
         $one_model = 1;
     }
      else {die "check parameter $ARGV[8] is not one???\n";}
 }


 my($file,$path_target,$path_seq,$target_output,$return_val,$path_read,$path_write);
 my(@tt,@files);

 opendir(DIR,"$dir_targets");
 @files=readdir(DIR);
 foreach $file (@files)
 {
	 if($file eq "." || $file eq "..")
	 {
		 next;
	 }
	 $target_output = $tem_output."/".$file;
	 -s $target_output || system("mkdir $target_output");
	 $path_target=$dir_targets."/".$file;
	 $path_seq = $dir_fasta."/".$file.".fasta";
         @tt = split(/\./,$file);
          $path_write= $dir_output."/".$tt[0].".ssp_similarity";
         if(check_it($path_write))
         {
              print "already processed $path_write\n";
              next;
         }

         if($one_model == 1)
         {
                 @tt = split(/\./,$file);
                 $path_seq = $dir_fasta."/".$tt[0].".fasta";
                  $return_val = system("perl $addr_feature_pl $path_target $path_seq $addr_LCS $dssp $dssp_pl $spx $target_output one");
                 if($return_val!=0) {die "perl $addr_feature_pl $path_target $path_seq $addr_LCS $dssp $dssp_pl $spx $target_output one fails!\n";}
                 $path_read = $target_output."/".$file.".tmp.ss_similarity";
                 $path_write= $dir_output."/".$tt[0].".ssp_similarity";
                  system("cp $path_read $path_write");
                 next;
         }



	 if(!-s $path_seq)
	 {
		 print "We don't find the sequence for target $file, check $path_seq\n";
		 next;
	 }
	 $return_val = system("perl $addr_feature_pl $path_target $path_seq $addr_LCS $dssp $dssp_pl $spx $target_output");
	 if($return_val!=0)
	 {
		 print "perl $addr_feature_pl $path_target $path_seq $addr_LCS $dssp $dssp_pl $spx $target_output fails!\n";
		 next;
	 }
     $path_read = $target_output."/".$file.".ss_similarity";
	 $path_write= $dir_output."/".$file.".ssp_similarity";
	 system("cp $path_read $path_write");
 }

 #system("rm -R $tem_output");
