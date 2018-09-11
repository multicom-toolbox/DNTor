#!/bin/sh
# Scripts/make_DeepProTa_prediction.sh  2mnjB  /exports/store1/jh7x3/DNTorsion/Inputs/2mnjB  DNN  /exports/store1/jh7x3/DNTorsion/output
echo ''
echo '  ###############################################################################'
echo '  #                                                                             #'
echo '  #     DeepProTa : Deep Learning Methods for Protein Torsion Angle Prediction  #'
echo '  #                                                                             #'
echo '  #     Copyright (C) 2016 -2024    Jie Hou, Haiou, and Jianlin Cheng           #'
echo '  #                                                                             #'
echo '  #     DeepProTa is free software: you can redistribute it and/or modify       #'
echo '  #     it under the terms of the GNU General Public License as published by    #'
echo '  #     the Free Software Foundation, either version 3 of the License, or       #'
echo '  #     (at your option) any later version.                                     #'
echo '  #                                                                             #'
echo '  #     DeepProTa is distributed in the hope that it will be useful,            #'
echo '  #     but WITHOUT ANY WARRANTY; without even the implied warranty of          #'
echo '  #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     See the        #'
echo '  #     GNU General Public License for more details.                            #'
echo '  #                                                                             #'
echo '  #     You should have received a copy of the GNU General Public License       #'
echo '  #     along with DeepProTa.     If not, see <http://www.gnu.org/licenses/>.   #'
echo '  #                                                                             #'
echo '  ###############################################################################'

echo ''

if [ "$#" -ne 4 ]; then
    echo "The number of parameters is not correct! Please run it as <./make_DeepProTa_prediction.sh sequence_name sequence_file(full path) modeltype(DNN,DRBM,DRNN,DReRBM) outputdir(full path)>"
    exit 1
fi
sequence_name=$1
sequence_file=$2
modeltype=$3
outputdir=$4

GLOBAL_PATH='/exports/store2/haiou/JieShare/Project_of_TorsionAngle/';

install_dir=${GLOBAL_PATH}
Scriptdir=${GLOBAL_PATH}/Scripts

cd ${Scriptdir}


# make prediction

if [ "$modeltype" == 'DNN' ];
then
   winsize='11'
elif [ "$modeltype" == 'DRBM' ];
then
   winsize='7'
elif [ "$modeltype" == 'DRNN' ];
then
   winsize='7'
elif [ "$modeltype" == 'DReRBM' ];
then
   winsize='3'
else
   echo "The model architecture is not satisfied! Only DNN, DRBM, DRNN, DReRBM is supported!"
   exit 1
fi


#use tools to generate features for sequence  
python ./get_Queryfile.py $sequence_name $sequence_file $install_dir $outputdir

#combine features
python ./combine_Features.py $sequence_name $outputdir $install_dir $sequence_file $modeltype


feature="${outputdir}/${sequence_name}.fea_${winsize}"
output="${outputdir}/${sequence_name}.${modeltype}.torsion"

echo "Setting feature to ${feature}"
echo "Setting output prediction to ${output}"
echo "Setting windowsize to ${winsize}"
python ./run_DeepProTa_prediction.py --feature  $feature --arch $modeltype  --output  $output
