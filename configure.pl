#!/usr/bin/perl -w
 use FileHandle; # use FileHandles instead of open(),close()
 use Cwd;
 use Cwd 'abs_path';

 #my($dir)=getcwd;
 #$dir=abs_path($dir);
 #print "We get this path:\n$dir\n";
 #exit(0);
#####################################################################################
#																					#
#configure.pl: to configure the installation of DeepQA								#
#																					#
#Author: Renzhi Cao                         										#
#Date: March, 2016																	#
#																					#
#####################################################################################

#####################################################################################
# Modified by: March , 2016															#
#																					#
#References:																		#
#																					#
#DeepQA: Master estimation of model accuracy with deep neural networks              #
#Renzhi Cao, Debswapna Bhattacharya, Jie Hou, Jianlin Cheng.		         		#
# 																					#
#Bioinformatics, submitted          												#
#####################################################################################


######################## !!! customize settings here !!! ############################
#																					#
# Set installation directory of DeepQA to your unzipped DeepQA directory            #
     
 $install_dir = "/your_path/DeepQA";
######################## !!! End of customize settings !!! ##########################

if($install_dir eq "/your_path/DeepQA")
{# user forgets to set the default path of DeepQA, try to solve this problem
    $install_dir = getcwd;
    $install_dir=abs_path($install_dir);
}





####################### Don't change the following code #############################


print "=============================================================================\n";
print "=                                                                           =\n";
print "= DeepQA : Master estimation of model accuracy with deep neural networks    =\n";
print "=                                                                           =\n";
print "=============================================================================\n";

if(!-s $install_dir)
{
	die "The DeepQA directory ($install_dir) is not existing, please revise the customize settings part inside the configure.pl, set the path as  your unzipped DeepQA directory\n";
}
if ( substr($install_dir, length($install_dir) - 1, 1) ne "/" )
{
        $install_dir .= "/";
}

print "checking whether the configuration file run in the installation folder ...";
$cur_dir = `pwd`;
chomp $cur_dir;
$configure_file = "$cur_dir/configure.pl";
if (! -f $configure_file || $install_dir ne "$cur_dir/")
{
        die "\nPlease check the installation directory setting and run the configure program in the installation directory of DeepQA.\n";
}
print " OK!\n";


$| = 1;
print "Checking the path of DeepQA ...";
# get the full path of script #
$script_path = $install_dir."script";
# get the full path of bin #
$bin_path = $install_dir."bin";

$tool_path=$install_dir."tools";
if(!-s $tool_path)
{
	die "Cannot find $tool_path, check the script setting of DeepQA folder : $install_dir\n";
}
if(!-s $script_path)
{
	die "Cannot find $script_path, check the script setting of DeepQA folder : $install_dir\n";
}
if(!-s $bin_path)
{
	die "Cannot find $bin_path, check the bin setting of DeepQA folder : $install_dir\n";
}

## check some scripts inside the script folder ##
$tmp = $script_path."/"."model2seq.pl";
#-s $tmp || die "Cannot find $tmp, check the setting of DeepQA folder\n";





print " Pass\n";



print "Generating DeepQA prediction file ...";
$DeepQA_sh = $bin_path."/"."DeepQA.sh";
$DeepQA_perl = $script_path."/"."DeepQA.pl";
open(SERVER_SH, ">$DeepQA_sh") || die "cannot write DeepQA shell script to $DeepQA_sh\n";
print SERVER_SH "#!/bin/sh\n#DeepQA prediction file for protein single-model quality assessment #\n";
print SERVER_SH "if [ \$# -ne 3 ]\n";
print SERVER_SH "then\n\techo \"need three parameters : path of fasta sequence, directory of input pdbs, directory of output\"\n\texit 1\nfi\n";
print SERVER_SH "perl $DeepQA_perl $script_path $tool_path \$1 \$2 \$3\n"; 
close SERVER_SH;
system("chmod 755 $DeepQA_sh");
system("chmod 755 $DeepQA_perl");
print " Done!\n";



#### change the permission of all scripts ####
opendir(DIR,"$script_path");
@all_scripts = readdir(DIR);
foreach $file (@all_scripts)
{
	if($file eq "." || $file eq "..")
	{
		next;
	}
	$path_1=$script_path."/".$file;
	system("chmod 755 $path_1");
}
##############################################
#### install pspro ####

print "Installing pspro2 ... ";
my($dir_pspro)=$tool_path."/"."pspro2/";


-s $dir_pspro || die "not exists $dir_pspro, please download the pspro2 (http://sysbio.rnet.missouri.edu/multicom_toolbox/tools.html) and put it at this path:$dir_pspro, or redownload DeepQA tool\n";
my($con_pspro) = $dir_pspro."/"."configure.pl";
my($DeepQA_pspro)=$dir_pspro."/"."DeepQA_configure.pl";
my($IN,$OUT,$line);
$OUT = new FileHandle ">$DeepQA_pspro";
$IN=new FileHandle "$con_pspro";
while(defined($line=<$IN>))
{
	chomp($line);
	if(substr($line,0,12) eq "\$install_dir")
	{
		print $OUT "\$install_dir = \"$dir_pspro\";\n"; 
	}
	elsif(substr($line,0,10) eq "\$nr_db_dir")
	{
		print $OUT "\$nr_db_dir = \"$tool_path/nr\";\n";
	}
	elsif(substr($line,0,11) eq "\$big_db_dir")
	{
		print $OUT "\$big_db_dir = \"$tool_path/big\";\n";
	}
	else
	{
		print $OUT $line."\n";
	}
}
$IN->close();
$OUT->close();


chdir("$dir_pspro");

system("chmod 755 $DeepQA_pspro");
$return_val=system("perl $DeepQA_pspro");
if($return_val)
{
	print "Installing tool pspro2 fails, check perl $DeepQA_pspro!\n";
	exit(0);
}

print "Done!\n";


print "Installing betacon ... ";
$dir_bete=$tool_path."/"."betacon/";
-s $dir_bete || die "not exists $dir_bete, please download the betacon(http://sysbio.rnet.missouri.edu/multicom_toolbox/tools.html) and put it at this path:$dir_bete, or redownload DeepQA tool\n";

my($con_beta) = $dir_bete."/"."configure.pl";
my($DeepQA_beta)=$dir_bete."/"."DeepQA_configure.pl";
$OUT = new FileHandle ">$DeepQA_beta";
$IN=new FileHandle "$con_beta";
while(defined($line=<$IN>))
{
	chomp($line);
	if(substr($line,0,12) eq "\$install_dir")
	{
		print $OUT "\$install_dir = \"$dir_bete\";\n"; 
	}
	elsif(substr($line,0,10) eq "\$pspro_dir")
	{
		print $OUT "\$pspro_dir = \"$dir_pspro/\";\n";
	}
	else
	{
		print $OUT $line."\n";
	}
}
$IN->close();
$OUT->close();


chdir("$dir_bete");
system("chmod 755 $DeepQA_beta");
$return_val=system("perl $DeepQA_beta");
if($return_val)
{
	print "Installing tool betacon fails, check perl $DeepQA_beta!\n";
	exit(0);
}

print "Done!\n";

print "Installing ModelEvaluator ... ";
$dir_modeleva=$tool_path."/"."model_eva1.0/";
chdir("$dir_modeleva");
-s $dir_modeleva || die "not exists $dir_modeleva, please download the ModelEvaluator (http://sysbio.rnet.missouri.edu/multicom_toolbox/tools.html) and put it at this path:$dir_bete, or redownload DeepQA tool\n";
my($con_modeleva) = $dir_modeleva."/"."configure.pl";
my($DeepQA_modeleva)=$dir_modeleva."/"."DeepQA_configure.pl";

$OUT = new FileHandle ">$DeepQA_modeleva";
$IN=new FileHandle "$con_modeleva";
while(defined($line=<$IN>))
{
	chomp($line);
	if(substr($line,0,12) eq "\$install_dir")
	{
		print $OUT "\$install_dir = \"$dir_modeleva\";\n"; 
	}
	else
	{
		print $OUT $line."\n";
	}
}
$IN->close();
$OUT->close();
chdir("$dir_modeleva");
system("chmod 755 $DeepQA_modeleva");
$return_val=system("perl $DeepQA_modeleva");
if($return_val)
{
	print "Installing tool ModelEvaluator fails, check perl $DeepQA_modeleva!\n";
	exit(0);
}

print "Done!\n";

print "Checking modeller ... ";
my($addr_mod913) = $tool_path."/"."modeller9.13/bin/mod9.13";
if (!-s $addr_mod913) {
	die "Please check $addr_mod913, you can download the modeller and install it by yourself if the current one in the tool folder is not working well, the key is MODELIRANJE.  please install it to the folder $tool_path/modeller9.13, with the file mod9.13 in the bin directory\n";
}

my($deep_mod913) = $tool_path."/"."modeller9.13/bin/modDeepQA";
$OUT = new FileHandle ">$deep_mod913";
$IN=new FileHandle "$addr_mod913";
while(defined($line=<$IN>))
{
        chomp($line);
        @ttt = split(/\=/,$line);

        if(@ttt>1 && $ttt[0] eq "MODINSTALL9v13")
        {
                print $OUT "MODINSTALL9v13=\"$tool_path/modeller9.13\"\n";
        }
        else
        {
                print $OUT $line."\n";
        }
}
$IN->close();
$OUT->close();
system("chmod 755 $deep_mod913");
my($modeller_conf) = $tool_path."/"."modeller9.13/modlib/modeller/config.py";
$OUT = new FileHandle ">$modeller_conf";
print $OUT "install_dir = r\'$tool_path/modeller9.13/\'\n";
print $OUT "license = \'MODELIRANJE\'";
$OUT->close();
system("chmod 755 $modeller_conf");
print "Done\n";

print "Installing scwrl4 ... ";
my($dir_scwrl) = $tool_path."/"."scwrl4";
-s $dir_scwrl || die "Cannot find $dir_scwrl, please re download the DeepQA tool or download scwrl4 by yourself\n";

my($addr_scwrl_orig) = $dir_scwrl."/"."Scwrl4.ini";
my($addr_scwrl_back) = $dir_scwrl."/"."Scwrl4.ini.back";
system("cp $addr_scwrl_orig $addr_scwrl_back");
my(@ttt);
$OUT = new FileHandle ">$addr_scwrl_orig";
$IN=new FileHandle "$addr_scwrl_back";
while(defined($line=<$IN>))
{
	chomp($line);
	@ttt = split(/\s+/,$line);
	
	if(@ttt>1 && $ttt[1] eq "FilePath")
	{
		print $OUT "\tFilePath\t=\t$dir_scwrl/bbDepRotLib.bin\n"; 
	}
	else
	{
		print $OUT $line."\n";
	}
}
$IN->close();
$OUT->close();
print "Done\n";

print "Installing spine_x .. ";
my($dir_spinex) = $tool_path."/"."spine_X";
-s $dir_spinex || die "Cannot find $dir_spinex, please re download the DeepQA tool or download spine_X by yourself\n";
my($addr_spinex_orig) = $dir_spinex."/"."spX.pl";
my($addr_spinx_new) = $dir_spinex."/"."DeepQA_spX.pl";

$OUT = new FileHandle ">$addr_spinx_new";
$IN=new FileHandle "$addr_spinex_orig";
while(defined($line=<$IN>))
{
	chomp($line);
	if(substr($line,0,14) eq "\$install_codir")
	{
		print $OUT "\$install_codir = \'$dir_spinex/code\';\n"; 
	}
	elsif(substr($line,0,14) eq "\$install_blast")
	{
		print $OUT "\$install_blast = \'$dir_spinex/blast-2.2.17\';\n"; 
	}
	elsif(substr($line,0,11) eq "\$install_nr")
	{
		print $OUT "\$install_nr = \"$tool_path/nr\";\n"; 
	}
	else
	{
		print $OUT $line."\n";
	}
}
$IN->close();
$OUT->close();
chdir("$dir_spinex");
system("chmod 755 $addr_spinx_new");
#$return_val=system("perl $addr_spinx_new");
#if($return_val)
#{
#	print "Installing tool spine_x fails, check perl $addr_spinx_new!\n";
#	exit(0);
#}

print "Done\n";


print "Install Qprob package ... ";
my($dir_qprob) = $tool_path."/"."qprob_package";
my($addr_qprob) = $dir_qprob."/"."configure.pl";
chdir("$dir_qprob");
system("chmod 755 $addr_qprob");
$return_val=system("perl $addr_qprob $tool_path");
if($return_val)
{
	print "Installing tool qprob fails, check perl $addr_qprob $tool_path!\n";
	exit(0);
}

print "Done\n";

print "Setting Deep Belief Network ...";
#my($addr_DBN_orig) = $tool_path."/"."DBN/run_prediction_orig.py";
#my($addr_DBN_new) = $tool_path."/"."DBN/Run_prediction_DeepQA.py";
#my($dir_DBN_lib) = $tool_path."/"."DBN/dn_libs/";
#if(!-s $addr_DBN_orig || !-s $dir_DBN_lib)
#{
#	die "Error, not found $dir_DBN_lib or $addr_DBN_orig, please double check if you download the full packages or remove some folders, you can download it again or contact the author of this package\n";
#}
#$OUT = new FileHandle ">$addr_DBN_new";
#$IN=new FileHandle "$addr_DBN_orig";
#while(defined($line=<$IN>))
#{
#	chomp($line);
#	if(substr($line,0,11) eq "GLOBAL_PATH")
#	{
#		print $OUT "GLOBAL_PATH=\'$dir_DBN_lib\';\n"; 
#	}
#	else
#	{
#		print $OUT $line."\n";
#	}
#}
#$IN->close();
#$OUT->close();
#system("chmod 755 $addr_DBN_new");
my($addr_DBN_new) = $tool_path."/"."DBN/run_DeepQA_prediction.pl";
system("chmod 755 $addr_DBN_new");
print "Done\n";

print "DeepQA successfully installed!\n\n";

print "DeepQA should be used as follows, please run:\n";
print "$bin_path/DeepQA.sh\n\n\n";


