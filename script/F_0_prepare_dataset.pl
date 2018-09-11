##########################################################################################################
#                  Function about creating 10 cross-validation train and test model for FANN             #
#																										 #
#										Renzhi Cao  													 #
#																										 #
#									    1/16/2013														 #
#																										 #
#																										 #
#									Revised at 2/4/2016 	                         					 #
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
 use List::Util qw/shuffle/;
  if (@ARGV != 4)
    { # @ARGV used in scalar context = number of args

	  print "\n\nThis script will add the target name at the end of each line, #T0515, ..., so later we may need this information!\n";

	  print("This program tries to create a whole training and testing dataset of QA for SVM. \n");
	  print "The input should be a folder, each subfold with the feature name has all targets predictions (T0515.score** ... ), the format is modelname prediction_score.... \n"; 
	  print("You should execute the perl program like this: perl $PROGRAM_NAME dir_all_models input_dir_all_features addr_SVM_with_comment addr_SVM\n");
	  print("\n********** example*********\n");
     print "perl $0 ../models ../normed_scores/predictions ../SVM_with_comment ../SVM\n";
		  exit(1) ;
    }
 my $starttime = localtime();
 print "\n The time started at : $starttime.\n";

 my($dir_models) = $ARGV[0];          # we make sure every model has a score
 my($features_dir)=$ARGV[1];               # the directory for all features 
 my($output)=$ARGV[2];                 # the output dataset
 my($SVM_output) = $ARGV[3];
#############################check the status#############################################################
 
 -e $features_dir || die "Not find input features $features_dir!\n";
 my($IN,$OUT,$OUT2,$key,$value,$line,$num_features,$i,$j,$path1,$path2,$real_path1,$real_path2,$size,$file,$target,$index_name,$interval,$k,$count,$count_2);
 my(@tem_split,@files,@targets,@name_list,@ranks,@shuffle_target,@tem2);
 my(%hash);
 my(%hash_feature);
 my($ii);
 #########################################################################################################
 #
 #
 #               1. Get the target name lists from score of filtered native and filtered prediction
 #
 #
 #########################################################################################################
 my(%hash_models) = ();
 opendir(DIR,"$dir_models");
 my(@models) = readdir(DIR);
 my($model_name);
 foreach $model_name (@models)
 {
     if($model_name eq "." || $model_name eq "..") {next;}
     $hash_models{$model_name} = 0;
 }

 #########################################################################################################
 #
 #
 #               2. Select the whole dataset
 #
 #
 #########################################################################################################
 

######## output #####################
     ########## output of training file###############
	 $path1=$output."_tmp";
	 if(-s $path1)
	 {
		 print "$path1 already exists!\n";
	 }
	 else
	 {
         open (File, "&gt;$path1");
         chmod (0777, $path1); 
         close (File);
	 }
	 $OUT=new FileHandle ">$path1";
	 defined($OUT) || die "Cannot open output file $path1 \n";
	 
         ########################################################################
         #    (2). Load the prediction score of each feature
         ########################################################################
		 $num_features=0;                    # the number of features
		 %hash_feature=();
		 opendir(DIR, "$features_dir");
         @files = readdir(DIR);
         foreach $file (@files)
         {
	        if($file eq '.' || $file eq '..')
	        {
	        	next;
	        }
			@tem_split=split(/\_/,$file);
			if(!exists $hash_feature{$tem_split[1]})
			{
				$hash_feature{$tem_split[1]}=$file;
			}
			$num_features++;
         }
=pod
print "now for real scores\n";
foreach $key (keys %hash) {print $key."\t".$hash{$key}."\n";
}
print "now for the features:\n";
foreach $key (keys %hash_feature) {print $key."\t".$hash_feature{$key}."\n";
}
=cut
                 %hash=();
		 for($ii=0;$ii<$num_features;$ii++)
		 {
			if(!exists $hash_feature{$ii})
			{
				print "The feature name in folder $features_dir should have the name like : 0_*,1_*,2_*...\n";
				exit(0);
			}		
			$real_path1=$features_dir."/".$hash_feature{$ii};
	                $file=$hash_feature{$ii};	
			opendir(DIR, "$real_path1");
			@targets=readdir(DIR);
			$real_path2="NULL";
			foreach $target (@targets) 
			{
				if($target eq '.' || $target eq '..')
				{
					next;
				}
                                @tem_split = split(/\./,$target);
                                if($tem_split[0] eq "DeepQA")
                                {
				    $real_path2=$real_path1."/".$target;
                                    last; 
                                }
			}
			if($real_path2 eq "NULL")
			{
				print "$name_list[$count] cannot find in $features_dir/$file\n";
				exit(0);
				next;
			}
#print "Check $real_path2\n";
                    foreach $model_name (keys %hash_models)
                    {
                        $hash_models{$model_name} = 0;
                    }

		    $IN=new FileHandle "$real_path2";
		    defined($IN) || die "Cannot open input file $real_path2";
            while ( defined($line = <$IN>))
            {
	          #read something here
	          chomp($line);
	          $line=~s/\s+$//;  # remove the windows character
              @tem_split=split(/\s+/,$line);
			  if(@tem_split<2)
			  {
				  next;
			  }
                          if(not exists $hash{$tem_split[0]})
                          {
                             $hash{$tem_split[0]} = -1;
                          }
			  $hash{$tem_split[0]}.="|".$tem_split[1];
                          if(exists $hash_models{$tem_split[0]}) {$hash_models{$tem_split[0]} = 1;}
		    }
                    # put the default 0.5 for all missing features #
                    foreach $model_name (keys %hash_models)
                    {
                         if($hash_models{$model_name}!=1)
                         {
                           if(not exists $hash{$model_name})
                           {
                             $hash{$model_name} = -1;
                           }  
                           $hash{$model_name}.="|"."0.5";                              
                         }
                    }
            $IN->close();
            } # end of for ii =0
		 ########################################################################
         #    (3). output real_value, prediction1, prediction2 ... into the tmp training file
         ########################################################################
		 foreach $key (keys %hash) 
		 {
			 @tem_split=split(/\|/,$hash{$key});

             if($hash{$key} =~ m/nan/)
			 {
				 print "We skip the nan value!\n";
				 next;
			 }

             if(@tem_split != $num_features+1)
			 {
				 print "\n".@tem_split.", now $num_features\nCheck here, $key , $hash{$key},  features miss: $real_path1, $real_path2 !\n";
				 #exit(0);
				 next;
			 }
             $size++;

			 for($value=0;$value<@tem_split;$value++)
			 {
				 $tem_split[$value]=sprintf("%.6f",$tem_split[$value]);
				 print $OUT $tem_split[$value]." ";
			 }
			 
			 print $OUT "#$name_list[$count] \t $key";          # add the comment

			 print $OUT "\n";
		 }
	 
     $OUT->close();
		 ########################################################################
         #    (4). Process the tmp training file to the real training file
         ########################################################################		 
	     $path2=$output;
	     if(-s $path2)
	     {
	    	 print "$path2 already exists!\n";
	     }
	     else
	     {
            open (File, "&gt;$path2");
            chmod (0777, $path2); 
            close (File);
    	 }
    	 $OUT=new FileHandle ">$path2";
         my($OUT2) = new FileHandle ">$SVM_output";
    	 defined($OUT) || die "Cannot open output file $path2 \n";

	#	 print $OUT $size."\t".$num_features."\t"."1"."\n";               # total number of data, number of features, number of output
         $key =0;
	     $IN= new FileHandle "$path1";
         defined($IN) || die "Cannot open input file $path1\n";
         while ( defined($line = <$IN>))
         {
	          #read something here
	          chomp($line);
	          $line=~s/\s+$//;  # remove the windows character
			  @tem_split=split(/\s+/,$line);
              print $OUT $tem_split[0]."\t";              # the real score
              print $OUT2 $tem_split[0]."\t";
			  $key = 0;
			  for($value=1;$value<@tem_split;$value++)
			  {
				  if(substr($tem_split[$value],0,1) eq "#")
				  {
					  $key = 1; 
				  }
				  if($key == 1)
				  {
					  print $OUT $tem_split[$value]."\t";
				  }
				  else
				  {
					  print $OUT $value.":".$tem_split[$value]."\t";
                                          print $OUT2 $value.":".$tem_split[$value]."\t";
				  }
			  }			  
			  print $OUT "\n";
                          print $OUT2 "\n";
		 }
		 $IN->close();

		 $OUT->close();
                 $OUT2->close();
   system("rm $path1");
