#! /usr/bin/perl -w
use FileHandle;

if(@ARGV<2)
{
    print "This script is going to collect the modelevaluator score of the predictions\n";
    print "For example:\n";
    print "perl $0 ../calculated_scores/feature_11_modeleva_predictions ../normed_scores/predictions/feature_11_modeleva_predictions\n";

    print "\nperl $0 ../Validation/1_calculated_scores/feature_11_modeleva_casp11_stage1 ../Validation/2_normed_scores/casp11_stage1/feature_11_modeleva_casp11_stage1_norm\n";
    print "perl $0 ../Validation/1_calculated_scores/feature_11_modeleva_casp11_stage2 ../Validation/2_normed_scores/casp11_stage2/feature_11_modeleva_casp11_stage2_norm\n";
    exit(0);
}

my($dir_input) = $ARGV[0];
my($dir_output) = $ARGV[1];

my(@files,@targets);
my($file,$path_model,$name,$path,$score,$model_name,$i);
my(@tem);
    $path = $dir_input;
    $path_out = $dir_output;
    $OUT= new FileHandle ">$path_out";
    opendir(DIR,$path);
    @targets = readdir(DIR);
    foreach $target (@targets)
    {
       if($target eq "." || $target eq "..")
       {next;}
       @tem = split(/\./,$target);
       if(@tem>2) 
       {
           if($tem[@tem-2] eq "dssp" || $tem[@tem-2] eq "set" || $tem[@tem-2] eq "inp") {print "Skip these not pdb models \n";next;}
       }
#print "check $target\n";
       if($tem[@tem-1] eq "model_eval_score")
       {
#print "in\n";
           $path_model = $path."/".$target;
           $score = "NULL";
#print "check $path_model\n";
           $IN=new FileHandle "$path_model";
           if(defined($line=<$IN>))
           {
               chomp($line);
               @tt = split(/\s+/,$line);
               if(@tt<3) {next;}
               $score = $tt[2];
           } 
           $IN->close();
#print "find score $score\n";
           $model_name = $tem[0];
           for($i=1;$i<@tem-1;$i++) {$model_name.=".".$tem[$i];}
           if($score eq "NULL") {$score=0.5;}
           print $OUT $model_name."\t".$score."\n";
       }     

    }
    $OUT->close();



