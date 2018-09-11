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


  if (@ARGV < 4)
    { # @ARGV used in scalar context = number of args
	  
          print "this script is going to calculate the qprob score\n";
	  print "\n************** Renzhi Cao *******************\n";
	  print "Input:\n";

	  print "0. Dir of all targets (each subfolder is a target tarball)!\n";
	  print "1. dir of fasta sequences \n";
	  print "2. address of qprob.pl in qprob package \n";
	  print "3. Dir of output\n";
          print "Optional, 4. one. for native pdbs\n"; 

	  print "\n**************** Renzhi Cao *****************\n";
	  print "\nFor example:\n";

print "\n********************** ab initio for validation ************************\n";
          print "perl $0 ../converted_deb_casp11 ../sequences /space2/rcrg4/transfer/tools_I_develop/qprob_package/qprob.pl ../1_calculated_scores/feature_15_qprob_DB\n";

	  exit(0);
	}

 
 my($dir_targets)=$ARGV[0];
 my($dir_fasta)=$ARGV[1];
 my($addr_feature_pl)=$ARGV[2];
 my($dir_output)=$ARGV[3];

 -s $dir_output || system("mkdir $dir_output");

 my($tem_output)=$dir_output."/TEM";
 -s $tem_output || system("mkdir $tem_output");

 my($one_model) = 0;
 if(@ARGV>4)
 {
     if($ARGV[4] eq "one")
     {
         print "setting one model ...\n";
         $one_model = 1;
     }
      else {die "check parameter $ARGV[3] is not one???\n";}
 }


 my($file,$path_target,$path_seq,$target_output,$return_val,$path_read,$path_write);
 my(@files,@tt);

 opendir(DIR,"$dir_targets");
 @files=readdir(DIR);
 foreach $file (@files)
 {
	 if($file eq "." || $file eq "..")
	 {
		 next;
	 }
	 $target_output = $tem_output."/".$file;
         
if(-s $target_output) {print "tem skip $target_output, later comment this one if needed\n";next;}

	 -s $target_output || system("mkdir $target_output");
	 $path_target=$dir_targets."/".$file;
	 $path_seq = $dir_fasta."/".$file.".fasta";
         @tt = split(/\./,$file);
         $path_write= $dir_output."/".$tt[0].".qprob";
         if(check_it($path_write))
         {
              print " $path_write already generated result\n";
              next;
         }
         if($one_model == 1)
         {
                 @tt = split(/\./,$file);
                 $path_seq = $dir_fasta."/".$tt[0].".fasta";
                 if(!-s $path_seq) {print "cannot find $path_seq";next;}
                  $return_val = system("perl $addr_feature_pl $path_seq $path_target $target_output one");
                 if($return_val!=0) {die "perl  $addr_feature_pl $path_seq $path_target $target_output one fails!\n";}
                 $path_read = $target_output."/".$file.".NOVEL_global_QA";
                 $path_write= $dir_output."/".$tt[0].".qprob";
                  system("cp $path_read $path_write");
                 next;
         }

	 if(!-s $path_seq)
	 {
		 print "We don't find the sequence for target $file, check $path_seq\n";
		 next;
	 }
	 $return_val = system("perl $addr_feature_pl $path_seq $path_target $target_output");
	 if($return_val!=0)
	 {
		 print "perl $addr_feature_pl $path_seq $path_target $target_output fails!\n";
		 next;
	 }
     $path_read = $target_output."/".$file.".NOVEL_global_QA";
	 $path_write= $dir_output."/".$file.".qprob";
	 system("cp $path_read $path_write");
 }

 #system("rm -R $tem_output");
