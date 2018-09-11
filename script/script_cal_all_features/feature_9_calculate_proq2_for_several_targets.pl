#! /usr/bin/perl -w
#
 require 5.003; # need this version of Perl or newer
 use English; # use English names, not cryptic ones
 use FileHandle; # use FileHandles instead of open(),close()
 use Carp; # get standard error / warning messages
 use strict; # force disciplined use of variables
 use Cwd;
 use Cwd 'abs_path';
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
 if(@ARGV<6)
 {
    print "For several targets, This script will use the Proq2 to do QA.\n";
    print "perl $0 addr_cal_proq2_one_target.pl addr_run_all_external.pl addr_score_rosetta addr_rosetta_database dir_targets dir_sequence dir_output\n";
    print "For example:\n";


print "\n********************** ab initio for validation ************************\n";
    print "perl $0 feature_9_calculate_proq2_for_one_target.pl /home/tool/rosetta_proq2/ProQ_scripts-master/bin/run_all_external.pl  /home/tool/rosetta_bin_linux_2015.39.58186_bundle/main/source/bin/score.linuxgccrelease /home/tool/rosetta_bin_linux_2015.39.58186_bundle/main/database ../converted_deb_casp11 ../sequences ../1_calculated_scores/feature_9_proQ2_DB\n";

    exit(0);
 }
 my($addr_cal_proq2)=$ARGV[0];
 my($addr_run_all)=$ARGV[1];
 my($score_rosetta)=$ARGV[2];
 my($rosetta_database)=$ARGV[3];
 my($dir_target)=$ARGV[4];
 my($dir_seq)=$ARGV[5];
 my($dir_out)=$ARGV[6];
 
 my($is_native) = 0;
 if(@ARGV>7)
 {
    $is_native = 1;
 } 
 -s $dir_out || system("mkdir $dir_out");

 my($file,$path_target,$name,$path_seq,$path_out);
 my(@files,@ttt);

 opendir(DIR,"$dir_target");
 @files = readdir(DIR);
 foreach $file (@files)
 {
    if($file eq "." || $file eq "..")
    {
       next;
    }
    $path_target = $dir_target."/".$file;
    $name = $file;
	$path_seq = $dir_seq."/".$name.".fasta";
    $path_out = $dir_out."/".$name.".proq2_score";
    if(-s $path_out)
    {
       print "exists $path_out. ..\n";
       next;
    }
    if($is_native == 1)
    {
        @ttt = split(/\./,$file);
        $path_seq = $dir_seq."/".$ttt[0].".fasta";
        $path_out = $dir_out."/".$ttt[0].".proq2_score";
        if(!-s $path_seq) { die "cannot find $path_seq\n";}
        if(check_it($path_out)) {next;}
        if(system("perl $addr_cal_proq2 $addr_run_all $score_rosetta $rosetta_database $path_target $path_seq $path_out one"))
        {
           print "perl $addr_cal_proq2 $addr_run_all $score_rosetta $rosetta_database $path_target $path_seq $path_out fails!\n";
     
        }
        next;
    }

	if(!-s $path_seq)
	{
		print "Not existing $path_seq!\n";
		next;
	}
    if(check_it($path_out)) {next;}
    if(system("perl $addr_cal_proq2 $addr_run_all $score_rosetta $rosetta_database $path_target $path_seq $path_out"))
    {
      print "perl $addr_cal_proq2 $addr_run_all $score_rosetta $rosetta_database $path_target $path_seq $path_out fails!\n";
      exit(0);
    }
 }




