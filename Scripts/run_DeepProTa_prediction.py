############################################################################
#
#     DeepProTa : Deep Learning Methods for Protein Torsion Angle Prediction
#
#     Copyright (C) 2016 -2024    Jie Hou, Haiou, and Jianlin Cheng
#
#     DeepProTa is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     DeepProTa is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with FUSION.     If not, see <http://www.gnu.org/licenses/>.
#
############################################################################


import sys
GLOBAL_PATH='/exports/store2/haiou/JieShare/Project_of_TorsionAngle/';
sys.path.insert(0, GLOBAL_PATH+'/Scripts/')


import numpy as np
from DNlibrary import loadmodel, importSVM, calc_hidden_probs,logistic_function, calc_linear_output,get_phi_psi
import os
import os.path
import optparse 


parser = optparse.OptionParser()
parser.add_option( '--feature', dest = 'feature',
    default = '',    # default empty!
    help = 'Feature data for protein sequence')
parser.add_option( '--arch', dest = 'arch',
    default = '',    # default empty!
    help = 'The architecture of deep learning to predict, options: DNN, DRBM, DRNN, DReRBM')
parser.add_option( '--output', dest = 'output',
    default = '',    # default empty!
    help = 'The output for torsion anlge prediction')
(options,args) = parser.parse_args()


DN_model = ""
if options.feature: # check feature file
    filepath = options.feature
    if not os.path.exists(filepath):
        print 'Error ! Couldn\'t find feature file. Exiting application...'
        sys.exit(1)
else:
    print 'Error ! Target feature is not defined correctly. Please try ...'
    print '    python2 Scripts/run_DeepProTa_prediction.py -h'
    print 'Exiting application...'
    sys.exit(1)

if options.arch: # check modeltype
    if options.arch == 'DNN':
        print "Setting architecture to ",options.arch
    elif options.arch == 'DRBM':
         print "Setting architecture to ",options.arch
    elif options.arch == 'DReRBM':
         print "Setting architecture to ",options.arch
    elif options.arch == 'DRNN':
         print "Setting architecture to ",options.arch
    else:
        print 'Error ! The architecture is not correct! Only support DNN,DRBM,DRNN,DReRBM. Exiting application...'
        sys.exit(1)
else:
    print 'Error ! Model architecture is not defined correctly. Please try ...'
    print '    python2 Scripts/run_DeepProTa_prediction.py -h'
    print 'Exiting application...'
    sys.exit(1)

if options.arch == 'DNN':
    DN_model = GLOBAL_PATH + '/Methods/Models/DNN/DBNmodel_0_util.dat'
elif options.arch == 'DRBM':
    DN_model = GLOBAL_PATH + '/Methods/Models/DRBM/DBNmodel_0_util.dat'
elif options.arch == 'DReRBM':
    DN_model = GLOBAL_PATH + '/Methods/Models/DReRBM/DBNmodel_0_util.dat'
elif options.arch == 'DRNN':
    DN_model = GLOBAL_PATH + '/Methods/Models/DRNN/DBNmodel_0_util.dat'
else:
    print 'Error ! The architecture is not correct! Only support DNN,DRBM,DRNN,DReRBM. Exiting application...'
    sys.exit(1)
    
data_file = options.feature #dataset
prediction = options.output

print ""
print "########################"
print "#  Model architecture: ", options.arch
print "#  Model file:         ", DN_model
print "#  Output prediction:  ", prediction
print "########################"
print ""


test_datafile = data_file
if not os.path.exists(test_datafile):
    print "File ",test_datafile, " doesn't exists\n"
test_data=importSVM(test_datafile);
 
test_feature = test_data
#print "test_label",test_label


if not os.path.exists(DN_model):
#if True:
   #dn.test_DeepLearningPro(test_data, target_data)
   print("Error: Model ",DN_model," doesn't exists!\n");
   exit(-1);

  
###############  load the models 
loadmodel(DN_model, globals())

outresult_test = prediction + '.tmp'
print('Running model on testing data <testing>: ')

data_l1 = test_feature
# quick sanity check
if (data_l1.shape[1] != l1_vh.shape[0]):  
   sys.stderr.write('There is a mismatch between weight dim and input dimension, please check if features or model architecture is set correctly!')
   sys.exit(1)

   

# DNN or DRBM
if options.arch == 'DNN' or options.arch == 'DRBM':

    data_l2 = calc_hidden_probs(data_l1, l1_vh, l1_hb)
    del data_l1

    data_l3 = calc_hidden_probs(data_l2, l2_vh, l2_hb)
    del data_l2

    data_l4 = calc_hidden_probs(data_l3, l3_vh, l3_hb)
    del data_l3

    final_score = calc_linear_output(data_l4, l4_vh, l4_hb)

    np.savetxt(outresult_test, final_score, fmt="%5f")

    get_phi_psi(outresult_test,prediction)

    os.remove(outresult_test)    
elif options.arch == 'DRNN' or options.arch == 'DReRBM':
    TrainingData=[]
    for nsample in range(data_l1.shape[0]): 
        X_time = data_l1[nsample]
        TrainingData.append(X_time.tolist())
    T = len(TrainingData)
    fea_dim=len(TrainingData[0])
    X_sub=np.zeros((T,fea_dim))
    for time_stamp in range(len(TrainingData)): 
        X_sub[time_stamp] = TrainingData[time_stamp]
    curr_probs = X_sub
    hidden_dim = l0_vh.shape[0]
    First_hidden_probs = np.zeros((T+1,hidden_dim))
    First_hidden_probs[-1] = np.zeros(hidden_dim)
    Finaloutput=[]
    for t in np.arange(T):
        prob = X_sub[t]
        First_hidden_probs[t] = logistic_function(np.dot(prob, l1_vh) + np.dot(First_hidden_probs[t-1], l0_vh) +l1_hb)
        data_l2 = First_hidden_probs[t]
        final_score = np.dot(data_l2, l2_vh)+l2_hb
        Finaloutput.append(final_score.tolist())
    
    # save it into a numpy array
    Finaloutput_array = np.vstack(Finaloutput)
    #print Finaloutput_array[0:5]
    np.savetxt(outresult_test, Finaloutput_array, fmt="%5f")

    get_phi_psi(outresult_test,prediction)

    os.remove(outresult_test)

print "Torsion angle prediction is finished!"
