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
 ## define some variables.
 my($author) = "Renzhi Cao";
 my($version) = "Version 1.0";
 my($reldate) = "May 23";
 if (@ARGV != 3)
    { # @ARGV used in scalar context = number of args
	  print("This program tries to use the tool scwrl to process the model, rebuild the side chain for each of them, this may improve the quality score of model evaluator tool! \n");
      print("You should execute the perl program like this: perl $PROGRAM_NAME input_target_dir addr_scwrl_exe output_dir\n");
	  print("For example:\n");
        print "perl $0 ../test/T0709 /space2/rcrg4/transfer/tools_I_develop/DeepQA/tools/scwrl4/Scwrl4 ../test/scwrl_pro_T0709\n";
        exit(1) ;
    }
 my($input_dir)=$ARGV[0];
 my($input_scwrl)=$ARGV[1];
 my($output_dir)=$ARGV[2];
 -e $input_dir || die "No input $input_dir\n";
 -e $input_scwrl || die "No input $input_scwrl\n";
 if(!-s $output_dir)
 {
	 system("mkdir $output_dir");
 }
 my($path_dir,$path1,$path2,$return_val,$path_input);
 ################ start to use scwrl to process each model ###############
opendir(DIR, "$input_dir");
my @files = readdir(DIR);
#foreach my $file (@files)
#{
#	if($file eq '.' || $file eq '..')
#	{
#		next;
#	}
##########do something to the file###################
	$path_dir=$output_dir;  # the output dir
	$path_input=$input_dir;
	opendir(DIR, "$path_input");
    my @targets = readdir(DIR);
    foreach my $target (@targets)
    {
	   if($target eq '.' || $target eq '..')
	   {
		   next;
	   }
	   $path1=$path_input."/".$target;
	   $path2=$path_dir."/".$target;
       $return_val=system("$input_scwrl -i $path1 -o $path2");
	   if($return_val!=0 || !-s $path2)
	   {
		   print "$input_scwrl -i $path1 -o $path2 fails!\n";
######### here maybe copy the old  to the new position??? #########
                  system("cp $path1 $path2");
#           sleep(10);
#		   exit(0);
	   }
    }
#}
