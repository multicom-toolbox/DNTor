#!/usr/bin/perl
######################################################################################################
# Modified by: March , 2016																			 #
#																								   	 #
#References:																						 #
#																									 #
#DeepQA: Master estimation of model accuracy with deep neural networks                               #
#																									 #
#																								     #
#Contact:                                                                                            #
#Renzhi Cao: rcrg4@mail.missouri.edu                                                                 #
#website: http://cactus.rnet.missouri.edu/DeepQA/                                                    #
#Bioinformatics, 2016    	      																     #
######################################################################################################

use Cwd;
use Cwd 'abs_path';
use FileHandle;



################### DO NOT change below ##################################
if (@ARGV < 5)
{
       print "perl $0 ../script/ ../tools ../test/T0709.fasta ../test/T0709 ../test/DeepQA_T0709\n";
	die "need input fasta file, target folder.\nExample:\nperl $0 ../script/ ../tools addr_fasta dir_target dir_output\nperl $0 ../script/ ../tools ../test/T0709.fasta ../test/T0709 ../test/DeepQA_T0709\nperl $0 ../script/ ../tools ../test/T0709.fasta ../test/T0709/server01_TS1 ../test/DeepQA_server01_only one";
}
$script_path = abs_path($ARGV[0]);
$tool_path = abs_path($ARGV[1]);

$addr_fasta = abs_path($ARGV[2]);
$dir_target = abs_path($ARGV[3]);
$dir_output = abs_path($ARGV[4]);

 my($one_model) = 0;
 my($tmp_targets) = $dir_output.".tmp";
 
 if(@ARGV>5)
 {
      if($ARGV[5] eq "one")
      {
          print "great, only acess one model\n";
          $one_model= 1;
		  -s $tmp_targets || system("mkdir $tmp_targets");
          system("cp $dir_target $tmp_targets/");
          $dir_target = $tmp_targets;
      }
 }


-s $dir_output || system("mkdir $dir_output");

$dir_output = abs_path($dir_output);
chdir($dir_output);
$addr_log = $dir_output."/"."LOGFILE";
open(data_log,">$addr_log");

$sequence = "NULL";
$IN = new FileHandle "$addr_fasta";
while(defined($line=<$IN>))
{
	chomp($line);
	if(substr($line,0,1) eq ">")
	{#head
		next;
	}
	if($sequence eq "NULL")
	{
		$sequence = $line;
	}
	else
	{
		$sequence.=$line;
	}
}
$IN->close();

print "Get sequence:\n$sequence\n";

chomp $sequence;

my($DeepQA_sequence) = $dir_output."/"."DeepQA.fasta";
system("cp $addr_fasta $DeepQA_sequence");

############   now use scwrl to process all models ############
#ab initio: generate model scores (for full length model only)   ####here is modelEvaluator
my($wen_path_scwrl)="$tool_path/scwrl4/Scwrl4";
my($betacon_dir) = "$tool_path/betacon";
my($eva_dir) = "$tool_path/model_eva1.0";

my($DeepQA_models)=$dir_output."/"."DeepQA";
-s $DeepQA_models || system("mkdir $DeepQA_models");

my($wen_path1,$wen_path2,$wen_path);
my(@wen_paths);
####### use scwrl to process the model #######
opendir(DIR, "$dir_target");
@wen_paths = readdir(DIR);
foreach $wen_path (@wen_paths)
{
        if($wen_path eq '.' || $wen_path eq '..')
        {
                next;
        }
        $wen_path1=$dir_target."/".$wen_path;
        $wen_path2=$DeepQA_models."/".$wen_path;
        $ren_return_val=system("$wen_path_scwrl -i $wen_path1 -o $wen_path2>/dev/null 2>&1");
        if($ren_return_val!=0)
        {
          print LOG "$wen_path_scwrl -i $wen_path1 -o $wen_path2 fails !\n";
        }
		if(!-s $wen_path2)
	    {# scwrl fails in this case, only CA atom is in the pdb, we directly copy it from original
			system("cp $wen_path1 $wen_path2");
		}
}

######## now you have $DeepQA_models and $DeepQA_sequence #########
my($script_process_dssp) = $script_path."/"."Feature_processing_parse_ss_use_dssp.pl";
my($addr_dssp) = $tool_path."/"."dsspcmbi";
my($addr_pro_dssp) = $script_path."/"."dssp2dataset.pl";
my($dir_DeepQA_dssp_pro) = $dir_output."/"."DeepQA_dssp_processed";
$ren_return_val = system("perl $script_process_dssp $DeepQA_models $addr_dssp $addr_pro_dssp $dir_DeepQA_dssp_pro>/dev/null 2>&1");
if($ren_return_val)
{
	die "Cannot process the model use dssp, please check : perl $script_process_dssp $DeepQA_models $addr_dssp $addr_pro_dssp $dir_DeepQA_dssp_pro\n";
}



####### generate each feature ###################
my($feature_out) = $dir_output."/"."f";                  # this is the Deep features
-s $feature_out || system("mkdir $feature_out");

print "Generating feature surface score ...";
my($feature1_out) = $feature_out."/"."feature_0_surface";
-s $feature1_out || system("mkdir $feature1_out");
my($feature1_output) = $feature1_out."/"."DeepQA.feature0_surface";
my($addr_script1) = $script_path."/"."feature_1_surface_score.pl";
$ren_return_val = system("perl $addr_script1 $dir_DeepQA_dssp_pro $feature1_output>/dev/null 2>&1 &");
if($ren_return_val)
{
	die "Cannot predict surface feature, please check : perl perl $addr_script1 $dir_DeepQA_dssp_pro $feature1_output\n";
}


print "\nGenerating feature Euclidean compact score ...";
my($feature8_out) = $feature_out."/"."feature_7_EC";
-s $feature8_out || system("mkdir $feature8_out");
my($feature8_output) = $feature8_out."/"."DeepQA.feature8_EC";
my($addr_script8) = $script_path."/"."feature_8_Euclidean_compact.pl";

$ren_return_val = system("perl $addr_script8 $DeepQA_models $feature8_output>/dev/null 2>&1");
if($ren_return_val)
{
	die "Cannot predict Euclidean Compact feature, please check : perl $addr_script7 $DeepQA_models $feature8_output \n";
}


print "\nGenerating feature OPUS score ...";
my($feature4_out) = $feature_out."/"."feature_3_OPUS";
-s $feature4_out || system("mkdir $feature4_out");
my($feature4_output) = $feature4_out."/"."DeepQA.feature4_OPUS";
my($addr_script4) = $script_path."/"."feature_4_OPUS_score.pl";
my($addr_OPUS) = $tool_path."/"."OPUS_PSP";

$ren_return_val = system("perl $addr_script4 $addr_OPUS $DeepQA_models $feature4_output>/dev/null 2>&1 &");
if($ren_return_val)
{
	die "Cannot predict OPUS feature, please check : perl $addr_script4 $addr_OPUS $DeepQA_models $feature4_output. Make sure the OPUS is installed\n";
}

print "\nGenerating feature RWplus score ...";
my($feature5_out) = $feature_out."/"."feature_4_RWplus";
-s $feature5_out || system("mkdir $feature5_out");
my($feature5_output) = $feature5_out."/"."DeepQA.feature5_RWplus";
my($addr_script5) = $script_path."/"."feature_5_RWplus_score.pl";
my($addr_RWplus) = $tool_path."/"."RWplus";

$ren_return_val = system("perl $addr_script5 $addr_RWplus $DeepQA_models $feature5_output>/dev/null 2>&1 &");
if($ren_return_val)
{
	die "Cannot predict RWplus feature, please check : perl $addr_script5 $addr_RWplus $DeepQA_models $feature5_output. Make sure the RWplus is installed\n";
}

print "\nGenerating feature Secondary structure penalty score ...";
my($feature7_out) = $feature_out."/"."feature_6_SP";
-s $feature7_out || system("mkdir $feature7_out");
my($feature7_output) = $feature7_out."/"."DeepQA.feature7_SP";
my($addr_script7) = $script_path."/"."feature_7_secondary_structure_penalty.pl";
my($addr_LCS) = $tool_path."/"."LCS";
my($addr_dssp) = $tool_path."/"."dsspcmbi";
my($addr_dssp_pro) = $script_path."/"."dssp2dataset.pl";
my($addr_spx) = $tool_path."/"."spine_X/DeepQA_spX.pl";

$ren_return_val = system("perl $addr_script7 $DeepQA_models $DeepQA_sequence $addr_LCS $addr_dssp $addr_dssp_pro $addr_spx $feature7_output>/dev/null 2>&1 &");
if($ren_return_val)
{
	die "Cannot predict Secondary structure penalty feature, please check : perl $addr_script7 $DeepQA_models $DeepQA_sequence $addr_LCS $addr_dssp $addr_dssp_pro $addr_spx $feature7_output. \n";
}


print "\nGenerating feature Qprob score ...";
my($feature9_out) = $feature_out."/"."f_8_Qp";
-s $feature9_out || system("mkdir $feature9_out");
my($feature9_output) = $feature9_out;
my($addr_script9) = $script_path."/"."feature_9_Qprob.pl";
my($addr_qprob) = $tool_path."/"."qprob_package/bin/Qprob.sh";

#print "perl $addr_script9 $DeepQA_models $DeepQA_sequence $addr_qprob $feature9_output\n";

$ren_return_val = system("perl $addr_script9 $DeepQA_models $DeepQA_sequence $addr_qprob $feature9_output>/dev/null 2>&1 &");
if($ren_return_val)
{
	die "Cannot predict Euclidean Compact feature, please check : perl $addr_script9 $DeepQA_models $DeepQA_sequence $addr_qprob $feature9_output \n";
}
print "\nGenerating feature GOAP score ...";
my($feature3_out) = $feature_out."/"."feature_2_GOAP";
-s $feature3_out || system("mkdir $feature3_out");
my($feature3_output) = $feature3_out."/"."DeepQA.feature3_GOAP";
my($addr_script3) = $script_path."/"."feature_3_GOAP_score.pl";
my($addr_goap) = $tool_path."/"."goap-alone/goap";
$ren_return_val = system("perl $addr_script3 $addr_goap $dir_target $feature3_output>/dev/null 2>&1 &");
if($ren_return_val)
{
	die "Cannot predict goap feature, please check : perl $addr_script3 $addr_goap $DeepQA_models $feature3_output. Make sure the goap is installed\n";
}

print "\nGenerating feature Dope score ...";
my($feature2_out) = $feature_out."/"."feature_1_dope";
-s $feature2_out || system("mkdir $feature2_out");
my($feature2_output) = $feature2_out."/"."DeepQA.feature2_dope_tmp";
my($addr_script2) = $script_path."/"."feature_2_dope_score.pl";
my($addr_modeller) = $tool_path."/"."modeller9.13/bin/modDeepQA";
$ren_return_val = system("perl $addr_script2 $addr_modeller $DeepQA_models $feature2_output>/dev/null 2>&1 &");
if($ren_return_val)
{
	die "Cannot predict dope feature, please check : perl $addr_script2 $DeepQA_models $feature2_output. Make sure you have installed modeller9.13 and set up the environment, check to run mod9.13, and see if you can find this command \n";
}

print "\nGenerating feature ModelEvaluator score ...";
my($feature6_out) = $feature_out."/"."feature_5_ModEva";
-s $feature6_out || system("mkdir $feature6_out");
my($feature6_output) = $feature6_out."/"."DeepQA.feature6_ModEva";
my($addr_script6) = $script_path."/"."feature_6_ModelEvaluator_score.pl";
my($addr_beta) = $tool_path."/"."betacon/bin/beta_contact_map.sh";
my($addr_modeva) = $tool_path."/"."model_eva1.0/bin/model_eval.sh";

$ren_return_val = system("perl $addr_script6 $DeepQA_models $DeepQA_sequence $addr_beta $addr_modeva $feature6_output>/dev/null 2>&1");
if($ren_return_val)
{
	die "Cannot predict Modevaluator feature, please check : perl $addr_script6 $DeepQA_models $addr_beta $addr_modeva $feature6_output. Make sure the ModelEvaluator is installed\n";
}


my($force_go) = $dir_output."/"."WAITING";
$OUT = new FileHandle ">$force_go";
print $OUT "Please delete me if you think some features may fail, but you are fine to use 0.5 as the default score, or you just don't want to wait\n";
$OUT->close();

my($check)=3;
my($tem) = $feature3_out."/"."DeepQA.feature3_GOAP.TMP";
while(1)
{

    if(!-s $tem && -s $feature3_output)
	{
		last;
	}
	sleep(rand(100));
	$check--;
	if($check>=0)
	{
	    print "waiting the goap generating scores - $feature3_output. You can delete the following file to continue next step, the missing feature values will be set to 0.5 as default, you may lose the accuracy: $force_go\n";
	}
	if(!-s $force_go)
    {
		print "not existing $force_go\n";
		last;
	}
}

my($feature4_OPUS) = $feature4_out."/"."DeepQA.feature4_OPUS";
my($feature8_Qprob) = $feature9_out."/"."DeepQA.Qprob_score";
my($tag)=0;
while(1)
{
	if(!-s $force_go)
    {
		last;
	}
	$tag=0;
	if(!-s $feature2_output)
	{
		print "2 waiting for $feature2_output ... | ";
		$tag=1;
	}
	if(!-s $feature4_OPUS)
	{
		print "4 waiting for $feature4_OPUS ... | ";
		$tag=1;
	}
        if(!-s $feature8_Qprob)
        {
                print "8 waiting for $feature8_Qprob ... | ";
                $tag=1;
        }
	if(!-s $feature5_output)
	{
		print "5 waiting for $feature5_output ... | ";
		$tag=1;
	}
	if($tag==0) 
	{
		last;
	}
	print "You can delete the following file to continue next step, the missing feature values will be set to 0.5 as default, you may lose the accuracy: $force_go\n";
    sleep(rand(30));
}

###### now normalize the scores #######

my($feature2_final) = $feature2_out."/"."DeepQA.feature2_dope";
my($norm_dope)=$script_path."/"."norm_feature_2_dope.pl";
system("perl $norm_dope $DeepQA_sequence $feature2_output $feature2_final");
system("rm $feature2_output");


my($feature4_OPUS_final) = $feature4_out."/"."DeepQA.feature4_OPUS_final";
my($norm_OPUS)=$script_path."/"."norm_feature_4_OPUS.pl";
system("perl $norm_OPUS $DeepQA_sequence $feature4_OPUS $feature4_OPUS_final");
system("rm $feature4_OPUS");

my($feature5_RW_final) = $feature5_out."/"."DeepQA.feature5_RWplus_final";
my($norm_RW)=$script_path."/"."norm_feature_5_RWplus.pl";
system("perl $norm_RW $DeepQA_sequence $feature5_output $feature5_RW_final");
system("rm $feature5_output");

my($feature7_final) = $feature7_output."/"."DeepQA.ss_similarity";
system("cp $feature7_final $feature7_out/");
system("rm -R $feature7_output");


my($feature6_final) = $feature6_out."/"."DeepQA.feature6_ModEva_final";
my($coll_modeva) = $script_path."/"."collect_feature_6_modeva.pl";
system("perl $coll_modeva $feature6_output $feature6_final");
system("rm -R $feature6_output");


my($feature3_final) = $feature3_out."/"."DeepQA.feature3_GOAP_final";
my($feature3_norm) = $script_path."/"."norm_feature_3_goap.pl";
system("perl $feature3_norm $DeepQA_sequence $feature3_output $feature3_final");
system("rm -R $feature3_output");



##### now convert the SVM format #######
my($SVM_com_out) = $dir_output."/"."SVM_with_comment";
my($SVM_out) = $dir_output."/"."SVM_format";
my($SVM_script) = $script_path."/"."F_0_prepare_dataset.pl";
my($SVM_log) = $dir_output."/"."F_0_prepare_dataset.log";

system("perl $SVM_script $dir_target $feature_out $SVM_com_out $SVM_out > $SVM_log");

####### now make predictions ########
print "\nNow make predictions by Deep Belief Network ...\n";
my($addr_predict) = $tool_path."/"."DBN/run_DeepQA_prediction.pl";
my($addr_DN_model) = $tool_path."/"."DBN/models";
my($addr_DN_out) = $dir_output."/"."Predictions.txt";
#my($python_path)="python2.7";  #revise this path if you have different version
$ren_return_val = system("perl $addr_predict $SVM_out $addr_DN_model $addr_DN_out");
if($ren_return_val)
{
    print "failed to run perl $addr_predict $SVM_out $addr_DN_model $addr_DN_out\n";
	exit(0);
}

my($addr_final_out)= $dir_output."/"."DeepQA_predictions.txt";
$IN = new FileHandle "$SVM_com_out";
@names=();
$i=0;
while(defined($line=<$IN>))
{
	chomp($line);
	@tem = split(/\s+/,$line);
	if(@tem>2)
	{
		$names[$i++]=$tem[@tem-1];
	}
}
$IN->close();
$count=0;
$IN = new FileHandle "$addr_DN_out";
$OUT = new FileHandle ">$addr_final_out";
while(defined($line=<$IN>))
{
	chomp($line);
	if($count>$i)
	{
		die "check the predictions $addr_DN_out and $SVM_com_out, the total number is not equal, should not happen\n";
	}
	if($line ne "")
	{
		print $OUT $names[$count++]."\t".$line."\n";
	}
}
$OUT->close();
$IN->close();

print "Done, the DeepQA has predicted the score at $addr_final_out\n";

# clean up #
system("rm -R $dir_output/DeepQA");
system("rm -R $dir_output/DeepQA.fasta*");
system("rm -R $dir_output/DeepQA_dssp_processed");
system("rm -R $dir_output/WAITING");
