import gzip
import cPickle as pickle
import time
import numpy as np
import sys
import os.path
import math

def loadmodel(fname, target_dict, verbose = True):
    fo = gzip.GzipFile(fname, 'rb')
    var_list = pickle.load(fo)
    if verbose:
        print "Load architecture: ", var_list
    for var in var_list:
        target_dict[var] = pickle.load(fo)
    fo.close()



def calc_hidden_probs(data, vh, hb):
   print "Calculating probs. of hidden layer, " + str(data.shape[0]) + " examples."
   probs = np.empty((data.shape[0], vh.shape[1]))
   #print "data.shape[0]: ",probs[data.shape[0]-1,:]

   print "Processing...",
   #for i in range(0,data.shape[0] ):
   #print "range(0,data.shape[0]): ",range(0,data.shape[0]-1)
   for i in range(0,data.shape[0]):
       probs[i,:] = 1. / (1 + np.exp(-(np.dot(data[i,:], vh) + hb)))
       if (i+1)%1000 == 0:
           print str(i+1) + "...",
           sys.stdout.flush()
   #print "data.shape[0]: ",probs[data.shape[0]-1,:]
   print "Done"
   return probs

 
def calc_linear_output(data, vh, hb):
   print "Calculating probs. of linear layer, " + str(data.shape[0]) + " examples."
   #print "Calculating probs. of linear layer, " + str(vh) + " examples."
   probs = np.empty((data.shape[0], vh.shape[1]))

   print "Processing...",
   for i in range(0,data.shape[0]):
       probs[i,:] =  np.dot(data[i,:], vh) + hb
       if (i+1)%1000 == 0:
           print str(i+1) + "...",
           sys.stdout.flush()

   print "Done"
   return probs


def importSVM(filename, delimiter=' ', comment='',skiprows=0, start=0, end = 0, dtype=np.float32):
    # Open a file
    file = open(filename, "r")
    print "Import file: ", file.name
    if skiprows !=0:
       dataset = file.read().splitlines()[skiprows:]
    if skiprows ==0 and start ==0 and end !=0:
       dataset = file.read().splitlines()[0:end]
    if skiprows ==0 and start !=0:
       dataset = file.read().splitlines()[start:]
    if skiprows ==0 and start !=0 and end !=0:
       dataset = file.read().splitlines()[start:end]
    else:
       dataset = file.read().splitlines()
    #print dataset
    newdata = []
    for i in range(0,len(dataset)):
        line = dataset[i]
        #print line
        if line[0] != comment:
           fea = line.split(delimiter)
           newline = []
           for j in range(0,len(fea)):
               if fea[j].find(':') >0 :
                   (num,val) = fea[j].split(':')
                   newline.append(float(val))
               
           #print newline    
           newdata.append(newline)
    #print newdata
    data = np.array(newdata, dtype=dtype)
    #print data
    print "Data transformation successed!\n"
    file.close()
    return data





def get_phi_psi(predict_result,output_angle):

    predict_file=open(predict_result,'r').readlines()
    output_file=open(output_angle,'w')
    
    output_file.write('PHI PSI\n')
    for i in xrange(len(predict_file)):
        #print i
        sin_phi=predict_file[i].strip().split()[0]
        cos_phi=predict_file[i].strip().split()[1]
        sin_psi=predict_file[i].strip().split()[2]
        cos_psi=predict_file[i].strip().split()[3]

        if (cos_phi==0):
            phi_value=90.0
        else:
            phi_value=float(np.arctan2(float(sin_phi),float(cos_phi)))*180/3.1415926
        if (cos_psi==0):
            psi_value=90.0
        else:
            psi_value=float(np.arctan2(float(sin_psi),float(cos_psi)))*180/3.1415926

        phi_value_2 = float("{0:.2f}".format(phi_value))
        psi_value_2 = float("{0:.2f}".format(psi_value))
        output_file.write(str(phi_value_2)+' '+str(psi_value_2)+'\n')    

    output_file.close()
  

def calculate_MAE(predicted_file,true_file,output_MAE):
    predict_result=open(predicted_file,'r').readlines()
    true_result=open(true_file,'r').readlines()
    output_file=open(output_MAE,'w')
    psi_sum=0.0
    phi_sum=0.0

    for x in xrange(1,len(true_result)):
        phi_P=float(predict_result[x].strip().split(' ')[0])
        psi_P=float(predict_result[x].strip().split(' ')[1])

        phi_E=float(true_result[x].strip().split(' ')[0])
        psi_E=float(true_result[x].strip().split(' ')[1])

        # perform evaluation
        psi_dev = psi_P - psi_E
        phi_dev = phi_P - phi_E
        if psi_dev < -180:
           psi_P = psi_P + 360
        if psi_dev > 180:
           psi_P = psi_P - 360
        if phi_dev < -180:
           phi_P = phi_P + 360
        if phi_dev > 180:
           phi_P = phi_P - 360
        psi_sum += abs(psi_P - psi_E)
        phi_sum += abs(phi_P - phi_E) 

    psi_MAE = psi_sum / len(true_result)
    phi_MAE = phi_sum / len(true_result)

    #np.savetxt(result_summary, out,  fmt='%.5f',delimiter=' ')
    output_file.write('*phi_MAE: '+str(phi_MAE)+' '+'*psi_MAE: '+str(psi_MAE)+'\n')  
    print "*phi_MAE is %.2f " % (phi_MAE)
    print "*psi_MAE is %.2f " % (psi_MAE)

def logistic_function(X):
    """
    Calculated objective function for 1 layer (ie, no hidden layer)
    discrimenator using cross entropy error.  W is a weight vector
    were the very last weight is taken to be the bias, X is a matrix
    of inputs for weights where each row is an example.  Targets is a
    1d array of 0s or 1s for class/target, should have same size as 
    number of examples (ie, X.shape[0])
    """
    Xdim = X.ndim
    if Xdim == 1 : # for logistic output, N*1 , input will be array in 1D, need transform it to 2D. 
       X = X.reshape(len(X),1)
    tmp1 = np.empty((X.shape[0], X.shape[1]))
    tmp1.fill(15.0)
    tmp2 = np.empty((X.shape[0], X.shape[1]))
    tmp2.fill(-15.0)
    #tmp1 = np.full((X.shape[0], X.shape[1]), 15.0)
    #tmp2 = np.full((X.shape[0], X.shape[1]), -15.0)
    X = np.where(~(X > 15), X,tmp1) #RuntimeWarning: overflow encountered in exp
    X = np.where(~(X < -15), X,tmp2)
    if Xdim == 1 :
       X = X.ravel()
    curr_probs = 1. / (1 + np.exp(-X))
    curr_probs = np.clip(curr_probs, 1e-15, (1-1e-15))

    return curr_probs
