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
 sub get_seq_from_atom;
 sub cal_distance($);
 sub get_coord($);
 
##############standard Amino Acids (3 letter <-> 1 letter)#######
my(%amino)=();
$amino{"ALA"} = 'A';
$amino{"CYS"} = 'C';
$amino{"ASP"} = 'D';
$amino{"GLU"} = 'E';
$amino{"PHE"} = 'F';
$amino{"GLY"} = 'G';
$amino{"HIS"} = 'H';
$amino{"ILE"} = 'I';
$amino{"LYS"} = 'K';
$amino{"LEU"} = 'L';
$amino{"MET"} = 'M';
$amino{"ASN"} = 'N';
$amino{"PRO"} = 'P';
$amino{"GLN"} = 'Q';
$amino{"ARG"} = 'R';
$amino{"SER"} = 'S';
$amino{"THR"} = 'T';
$amino{"VAL"} = 'V';
$amino{"TRP"} = 'W';
$amino{"TYR"} = 'Y';
###################################################################
  
  if (@ARGV != 2)
    { # @ARGV used in scalar context = number of args
	  print "This feature comes from paper <<capturing native/native like structures with a physico-chemical metric in protein folding>>, 2013. biochim,biophys.\n";
          print "Revised by 2/19/2014, use 1 - score, so the correlation is positive!\n";
	  print"This script process all protein target, each subfolder in a input directory is a protein target, we evaluate the quality of model based on the pairwise summation of all euclidean distance between atoms!\n";
	  print "\n************** Renzhi Cao *******************\n";
	  print "Input:\n";
	  print "0. Dir of input target folder, for make prediction!\n";
	  print "1. Directory of output. \n";
      

	  print "\n**************** Renzhi Cao *****************\n";
	  print "\nFor example:\n";
print "\n********************** ab initio for validation ************************\n";
          print "perl $0 ../converted_deb_casp11 ../1_calculated_scores/feature_14_euclid_DB\n";


	  exit(0);
	}
############### the copy right for each script #############################
print  "REM       \n";
print  "REM          *********************************************************************************\n";
print  "REM          *                                                                               *\n";
print  "REM          *                                                                               *\n";
print  "REM          *      Developed By :   Renzhi Cao (rcrg4\@mail.missouri.edu)                   *\n";
print  "REM          *      Copyright    :   Dr. Jianlin Cheng's BDML Laboratory                     *\n";
print  "REM          *      Release Date :   June 13, 2013                                          *\n";
print  "REM          *      Vesion       :   1.0                                                     *\n";
print  "REM          *                                                                               *\n";
print  "REM          *********************************************************************************\n";
print  "REM       \n";
############################################################################

 my($addr_input)=$ARGV[0];

 my($dir_output)=$ARGV[1];

 -s $addr_input || die "cannot open input $addr_input\n";

 -s $dir_output || system("mkdir $dir_output");

 my($IN,$OUT,$line,$name,$path_out,$path_target,$target,$file,$path_pdb,$return_val,$len,$ideal_radius,$pair_wise_distance,$quality);
 my(@tem_split,@files,@targets);
 $| = 1;
 opendir(DIR,"$addr_input");
 @files=readdir(DIR);
 foreach $file (@files)
 {

	if($file eq '.' || $file eq '..')
	{
	  next;
    }
	$path_target=$addr_input."/".$file;          # the target folder
print "Processing $path_target ...\n";
	$name=$file;
    $path_out = $dir_output."/".$file.".pairwise_distance";      # the output file

#if(-s $path_out)
#{
#	next;
#}

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
		 $path_pdb=$path_target."/".$target;
		 #$len=get_seq_from_atom($path_pdb);        # get the sequence length
		 #$ideal_radius=2.2 * ($len**0.38);           # this is the ideal radius of gyration
         $pair_wise_distance=cal_distance($path_pdb); # calculate the pairwise distance of all atoms from the pdb
print "Model : $path_pdb \n";

#print "Real radius gyration is $actual_radius, ideal radius of gyration is $ideal_radius.\n";

         $quality=$pair_wise_distance;
	 
         $quality=1-$quality;

         print $OUT $target."\t".$quality."\n";

	}
    $OUT->close();
	
 }


 #read sequence length from atom file
sub get_seq_from_atom
{
	#assume the atom file exists
	my $file = $_[0];

	open(ATOM, $file) || die "can't read atom file: $file\n";
	my @atoms = <ATOM>;
	close ATOM; 
	my $prev = -1;
	my $seq = ""; 
	while (@atoms)
	{
		my $text = shift @atoms;
		if ($text =~ /^ATOM/)
		{
			#get aa name
			#get position
			my $res = substr($text, 17, 3); 
			$res = uc($res); 
			$res =~ s/\s+//g;
			my $pos = substr($text, 22, 4);

			#if 3-letter, convert it to 1-letter.
			if (length($res) == 3)
			{
				if (exists($amino{$res}) )
				{
					$res = $amino{$res}; 
				}
				else
				{
					$res = "X"; 
					print "$file: resudie is unknown, shouldn't happen.\n"; 
				}
			}
			if ($pos != $prev)
			{
				$seq .= $res; 
				$prev = $pos; 
			}
		}
	}
	return length($seq); 
}

sub cal_distance($)
{# calculate the pairwise distance between all atoms
	my($path_pdb)=@_;
	my($IN,$line,$j);
	my(@cor_X)=();
	my(@cor_Y)=();
	my(@cor_Z)=();
	my($index)=0;
      	my(@tem_cor)=();
        my(@aa)=();              # the amino acid for the input pdb
        my($res);


	my($i)=0;
	#### read pdb and calculate the radius of gyration ####
	 $IN=new FileHandle "$path_pdb";
	 while(defined($line=<$IN>))
	 {
		 chomp($line);
		 $line=~s/\s+$//;
	     if(substr($line,0,6) ne "ATOM  ")
		 {# this is not ATOM line
			 next;
		 }
         ###### check whether this is CA #########
		 if(substr($line,13,2) ne "CA" && substr($line,12,2) ne "CA")
		 {# the reason I put substr($line,12,2) is because some pdbs, the CA is not at the proper location.
			 next;
		 }
                 $res = substr($line,17,3);
                 $res=uc($res);
                 $res=~s/\s+//g;
                 if(not exists $amino{$res})
                 {
                     print "not exists residue $res, check the pdb!\n";
                     next;
                 }
                 $aa[$index]=$amino{$res};       # convert 3 to 1 amino acid
		 @tem_cor=get_coord($line);      # get the coordinate
		 $cor_X[$index]=$tem_cor[0];
		 $cor_Y[$index]=$tem_cor[1];
		 $cor_Z[$index]=$tem_cor[2];
		 $index++;

	 }
	 $IN->close();
	 if($index<=2)
	 {
		 print "There is no atoms in input pdb $path_pdb, check here !!!!\n";
		 return 0;
		 #exit(0);
		 #return 1;
	 }

     ###### calculate the pairwise distance between all unique pair atoms ######
         my(%pair_dis)=();               # the key is the unique amino acid pair, could be AA,AB, ... . The value is total distance and the frequency, finally get the average for each amino acid pair
         my(%straight_pair_dis)=();       # this is for calculating the dis for straight protein of the same sequence

         my($key,$value,$for_str);
         

	 my($radius)=0;
	 for($i=0;$i<$index-1;$i++)
	 {
		 for($j=$i+1;$j<$index;$j++)
		 {
                    $key=$aa[$i].$aa[$j];    # the key
		    $radius= sqrt( ($cor_X[$i]-$cor_X[$j])*($cor_X[$i]-$cor_X[$j]) + ($cor_Y[$i]-$cor_Y[$j])*($cor_Y[$i]-$cor_Y[$j]) + ($cor_Z[$i]-$cor_Z[$j])*($cor_Z[$i]-$cor_Z[$j]) );   # the euclidean distance for the aa pair
                    if(not exists $pair_dis{$key})
                    { # this is the first time we see this unique key
                        $pair_dis{$key} = $radius."|"."1"; 
                        $for_str=($j-$i)*3.8;
                        $straight_pair_dis{$key} = $for_str."|"."1";
                    }
                    else
                    {# already exists
                        @tem_split = split(/\|/,$pair_dis{$key});
                        $tem_split[0]+=$radius;
                        $tem_split[1]++; 
                        $value=$tem_split[0]."|".$tem_split[1]; 
                        $pair_dis{$key}=$value;        # update the value
                      
                        @tem_split = split(/\|/,$straight_pair_dis{$key}); 
                        $for_str=($j-$i)*3.8;
                        $tem_split[0]+=$for_str; 
                        $tem_split[1]++;
                        $for_str=$tem_split[0]."|".$tem_split[1];
                        $straight_pair_dis{$key} = $for_str;
                    }
		 }
	 }
     ######## calculate the avarage euclidean distance for each unique aa pair , and take the summation##########
     $radius=0;
     $for_str=0;      # this is for straight pdb
     foreach $key (keys %pair_dis)
     {
         @tem_split=split(/\|/,$pair_dis{$key});
         if(@tem_split<2)
         {
            next;
         }
         $value=$tem_split[0]/$tem_split[1];       # get the avarage value
         
         $radius+=$value;                  # get the summation of average

         @tem_split=split(/\|/,$straight_pair_dis{$key});
         $for_str+=$tem_split[0]/$tem_split[1];
 
     }
     if($for_str == 0)
     {
         return 0;
     }     
     $radius/=$for_str;
     

     return $radius;

}

 sub get_coord($)
 {
	 my($line)=@_;
	 my(@coor)=();
	 my($num);
	 $num=substr($line,30,8);
	 $coor[0]=sprintf("%.3f",$num);
	 $num=substr($line,38,8);
	 $coor[1]=sprintf("%.3f",$num);
	 $num=substr($line,46,8);
	 $coor[2]=sprintf("%.3f",$num);

	 return @coor;
     
 }

