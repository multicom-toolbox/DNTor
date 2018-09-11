##########################################################################################################
#                  Function about calculating the secondary structure for all pdb          				 #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    10/4/2012														 #
#																										 #
#																										 #
#									Revised at 10/4/2012	                         					 #
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
	  print("This script tries to use DSSP to parse the secondary structure for all models, the input folder has folders for each target, and the output folder will have the sequence and the secondary structure\n");
	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir_input_pdbs address_DSSP addr_dssp2dataset.pl dir_output\n");
      print("\n********** example******\n");
          print "perl $0 ../converted_deb_casp11 /space2/rcrg4/transfer/tool/dsspcmbi dssp2dataset.pl ../dssp_parsed_DB\n";
  	  exit(1) ;
    }
 my $starttime = localtime();
 print "\n The time started at : $starttime.\n";
 my($input_dir)=$ARGV[0];
 my($addr_dssp)=$ARGV[1];
 my($addr_dssp_parser)=$ARGV[2];
 my($output_dir)=$ARGV[3];
 -e $input_dir || die "$input_dir not exists!\n";
 -e $addr_dssp || die "$addr_dssp not exists!\n";
 -e $addr_dssp_parser || die "$addr_dssp_parser not exists!\n";

######################################################################################################## #\
 my($IN,$OUT,$i,$line,$read_folder,$target,$write_dssp,$return_val,$read_model,$dssp_parsed);
 my(@tem_split,@targets,@proq_scores,@names,@files);

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

print "Input model folder is $input_dir!\n";
print "Check succeed, temporary file created! Start ...\n";

my($target);
my(@targets);
my($path_write,$path_target);


	$path_target=$input_dir;
	$path_write=$output_dir;
	-s $path_write || system("mkdir $path_write");

	opendir(DIR, "$path_target");
	@files=readdir(DIR);
	foreach my $file (@files)
	{
		if($file eq '.' || $file eq '..')
		{
			next;
		}
		

	    $read_folder=$path_target."/".$file;
        
		
	    $write_dssp=$path_write."/".$file.".dsspout";         # the path for store the dssp output of this target $file
		if(-s $write_dssp)
		{# already processed
			print "$file already processed, since $write_dssp exists! next ...\n";
			next;
		}
		else
		{
		     open (File, "&gt;$write_dssp");
			 chmod (0777, $write_dssp);
	         close (File);		
		}
	    $return_val=system("$addr_dssp $read_folder $write_dssp");
		if($return_val!=0)
		{
			print "$addr_dssp $read_folder $write_dssp fails!\n";
#		exit(0);
			next;
		}
		$dssp_parsed=$path_write."/".$file.".dssp_parsed";  # the path for parsed dssp output result
		if(-s $dssp_parsed)
		{# already processed
			print "$file already processed, since $dssp_parsed exists! next ...\n";
			next;
		}
		else
		{
			 open (File, "&gt;$dssp_parsed");
		     chmod (0777, $dssp_parsed);
		     close (File);		
		}
		$return_val=system("perl $addr_dssp_parser $write_dssp $dssp_parsed");
		if($return_val!=0)
		{
			print "perl $addr_dssp_parser $write_dssp $dssp_parsed fails!\n";
			exit(0);
		}
		system("rm $write_dssp");           # remove the original dssp output, we can just use the parsed dssp output .
	}#end inside for each
 my $endtime = localtime();
 print  "\nThe time ended at : $endtime.\n";
