import sys
import os
import os.path
import math
import shutil
import numpy 
import time

#global winsize
winsize=3
#global half_winsize
half_winsize=int(winsize/2)

def set_windowsize_by_modeltype(type):
    global winsize
    global half_winsize
    print "\n###############################################################################"
    if modeltype == 'DRNN':
        winsize = 7 
        half_winsize=int(winsize/2) 
        print "# Setting architecture as DRNN (Deep Recurrent Neural Network), window size is set to 7!" 
    elif modeltype == 'DRBM':
        winsize = 7
        half_winsize=int(winsize/2)
        print "# Setting architecture as DRBM (Deep Restricted Boltzmann Machine), window size is set to 7!" 
    elif modeltype == 'DReRBM':
        winsize = 3
        half_winsize=int(winsize/2)
        print "# Setting architecture as DReRBM (Deep Recurrent Restricted Boltzmann Machine), window size is set to 3!" 
    elif modeltype == 'DNN':
        winsize = 11
        half_winsize=int(winsize/2)
        print "# Setting architecture as DNN (Deep Neural Network), window size is set to 11!" 
    else:
        print "DeepProTa only supports winsize to 3(for DReRBM)or 7(for DRBM and DRNN) or 11 (for DNN)"
        exit(-1)

def window(feature):
    # apply the windowing to the input feature
    feat = numpy.array(feature)
    #print 'winsize='+str(half_winsize)

    output = numpy.concatenate([numpy.vstack([feat[0]]*half_winsize), feat])

    output = numpy.concatenate([output, numpy.vstack([feat[-1]]*half_winsize)])
    output = [numpy.ndarray.flatten(output[i:i+2*half_winsize+1]).T for i in range(0,feat.shape[0])]
    return output

def window_data(*feature_types):
    n = len(feature_types[0])
    features = numpy.empty([n,0])

    for feature_type in feature_types:
        test = numpy.array(window(feature_type))
        features = numpy.concatenate((features, test), axis=1)

    return features
#======================================================================================================================

#========================================get PSSM======================================================================
def get_input_PSSM(pssm_file):
    idx_res = (0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18)
    fp = open(pssm_file, 'r')
    lines = fp.readlines()
    fp.close()

    pssm = []
    # iterate over the pssm file and get the needed information out
    for line in lines:
        split_line = line.split()
        if (len(split_line) == 22) and (split_line[0] != '#'):
            pssm_temp = [float(i) for i in split_line[2:22]]
            pssm.append([pssm_temp[k] for k in idx_res])

    return pssm

def normalize_PSSM(input_PSSM,final_output,min_max_dir,pdb_name,max_pssm,min_pssm):
    input_dir=input_PSSM
    output_dir=final_output

    max_list=open(max_pssm,'r').readlines()
    min_list=open(min_pssm,'r').readlines()

    feature_list=open(input_dir,'r').readlines()
    output_feature_file=open(output_dir+'/'+pdb_name+'_pssm_own_new','w')

    #print len(feature_list)
    #print output_dir+'/'+pdb_name+'_pssm_own_new'
    for i in xrange(len(feature_list)):
        for j in xrange(len(feature_list[i].split(' '))-1):
            feature_value=float(float(feature_list[i].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))/(float(max_list[0].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))
            if(float(feature_value)>1):
                feature_value=1.0
            if(float(feature_value)<0.0):
                feature_value=0.0
            output_feature_file.write(str(feature_value)+' ')
        output_feature_file.write('\n')

def get_PSSM(feature_number,normalize,in_dir,out_dir,normalize_dir,final_output_dir,winsize,sequence_name):

    if_normalize=normalize
    input_dir=in_dir
    output_dir=out_dir

    max_pssm=normalize_dir+'/max_pssm'
    min_pssm=normalize_dir+'/min_pssm'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    pdb_name=sequence_name
    #print pdb_name
    pssm_file=input_dir
    output_file=open(output_dir+'/'+pdb_name+'.PSSM_own_new','w')
    pssm = get_input_PSSM(pssm_file)
    #print aa
    input_feature = window_data(pssm)

    for i in xrange(len(input_feature)):
        for j in xrange(len(input_feature[i])):
             output_file.write(str(input_feature[i][j]) + ' ')
        output_file.write('\n')
    output_file.close()
    #f1.close()
    #normalize PSSM
    normalize_PSSM(output_dir+'/'+pdb_name+'.PSSM_own_new',final_output_dir,out_dir,pdb_name,max_pssm,min_pssm)

#=========================================finish PSSM=================================================================

#==========================================get ASA====================================================================

def get_input_asa(asa_file,feature_number):
    fp = open(asa_file, 'r')
    lines = fp.readlines()
    fp.close()

    asa=[]

    for line in lines:
        split_line=line.split()

        if(len(split_line) == int(feature_number)):
            asa_tmp= [int(i) for i in split_line[0:int(feature_number)]]
            asa.append(asa_tmp)

    return asa

def get_ASA(feature_number,normalize,in_dir,out_dir,final_output_dir,winsize,sequence_name):

    if_normalize=normalize
    input_dir=in_dir
    output_dir=out_dir

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    fasta_name=sequence_name
    #=encode solvent accessibility

    input_feature=open(input_dir).readlines()
    output_feature=open(output_dir+'/'+fasta_name+'_encode.scratchasa','w')

    for i in xrange(len(input_feature[1])):
        if (input_feature[1][i] == 'e'):
            output_feature.write('1\n')

        if (input_feature[1][i] == '-'):
            output_feature.write('0\n')    

    output_feature.close()


    #window solvent accessibility

    pdb_name=sequence_name
    #print pdb_name
    output_file=open(final_output_dir+'/'+pdb_name+'_'+str(feature_number)+'_asa','w')

    asa_file=output_dir+'/'+pdb_name+'_encode.scratchasa'

    asa=get_input_asa(asa_file,feature_number)
    input_feature_final = window_data(asa)

    for i in xrange(len(input_feature_final)):
        for j in xrange(len(input_feature_final[i])):
             output_file.write(str(input_feature_final[i][j]) + ' ')
        output_file.write('\n')

    output_file.close()

#============================================finish ASA===============================================================

#==========================================get ss=====================================================================
def get_input_ss(ss_file,feature_number):

    #print ss_file
    fp = open(ss_file, 'r')
    lines = fp.readlines()
    fp.close()

    ss=[]
    #print "lines number="+str(len(lines))
    for line in lines:
        split_line=line.split()

        if(len(split_line) == int(feature_number)):
            ss_tmp= [(i) for i in split_line[0:int(feature_number)]]
            ss.append(ss_tmp)
    return ss

def get_SS_8class(feature_number,normalize,in_dir,out_dir,final_output_dir,winsize,sequence_name):

    if_normalize=normalize
    input_dir=in_dir
    output_dir=out_dir

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    #encode secondary structure
    fasta_name=sequence_name
    #print fasta_name
    input_feature=open(input_dir,'r').readlines()
    output_feature=open(output_dir+'/'+fasta_name+'_encode.scratchss','w')

    for i in xrange(len(input_feature[1])):
        if (input_feature[1][i] == 'G'):
            output_feature.write('1 0 0 0 0 0 0 0\n')

        if (input_feature[1][i] == 'H'):
            output_feature.write('0 1 0 0 0 0 0 0\n')    

        if (input_feature[1][i] == 'I'):
            output_feature.write('0 0 1 0 0 0 0 0\n')    

        if (input_feature[1][i] == 'T'):
            output_feature.write('0 0 0 1 0 0 0 0\n')

        if (input_feature[1][i] == 'E'):
            output_feature.write('0 0 0 0 1 0 0 0\n')    

        if (input_feature[1][i] == 'B'):
            output_feature.write('0 0 0 0 0 1 0 0\n')    

        if (input_feature[1][i] == 'S'):
            output_feature.write('0 0 0 0 0 0 1 0\n')

        if (input_feature[1][i] == 'C'):
            output_feature.write('0 0 0 0 0 0 0 1\n')    

    output_feature.close()    


    #window for secondary structure  

    pdb_name=sequence_name
    #print pdb_name
    output_file=open(final_output_dir+'/'+pdb_name+'_'+str(feature_number)+'_ss','w')

    ss_file=output_dir+'/'+pdb_name+'_encode.scratchss'
    ss=get_input_ss(ss_file,feature_number)
    input_feature_final = window_data(ss)

    for i in xrange(len(input_feature_final)):
        for j in xrange(len(input_feature_final[i])):
             output_file.write(str(input_feature_final[i][j]) + ' ')
        output_file.write('\n')

    output_file.close()
#===========================================finish ss=================================================================


#==================================get 7PC============================================================================
def normalize_7PC(input_7PP,final_output,pdb_name,max_7PC,min_7PC):
    input_dir=input_7PP
    output_dir=final_output

    max_list=open(max_7PC,'r').readlines()
    min_list=open(min_7PC,'r').readlines()

    feature_list=open(input_dir,'r').readlines()
    output_feature_file=open(output_dir+'/'+pdb_name+'_7pc','w')
    for i in xrange(len(feature_list)):
        for j in xrange(len(feature_list[i].split(' '))-1):
            feature_value=float(float(feature_list[i].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))/(float(max_list[0].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))
            if(float(feature_value)>1):
                feature_value=1.0
            if(float(feature_value)<0.0):
                feature_value=0.0
            output_feature_file.write(str(feature_value)+' ')
        output_feature_file.write('\n')


def get_7PC(feature_number,normalize,in_dir,out_dir,final_output_dir,winsize,normalize_dir,sequence_name):

    if_normalize=normalize
    input_dir=in_dir
    output_dir=out_dir

    max_7PC=normalize_dir+'/max_7PC'
    min_7PC=normalize_dir+'/min_7PC'

    # define the dictionary with the phys properties for each AA
    phys_dic = {            'A': [1.28, 0.05, 1.00, 0.31, 6.11, 0.42, 0.23],
                            'C': [1.77, 0.13, 2.43, 1.54, 6.35, 0.17, 0.41],
                            'D': [1.60, 0.11, 2.78, -0.77, 2.95, 0.25, 0.20],
                            'E': [1.56, 0.15, 3.78, -0.64, 3.09, 0.42, 0.21],
                            'F': [2.94, 0.29, 5.89, 1.79, 5.67,0.30, 0.38],
                            'G': [0.00, 0.00, 0.00, 0.00, 6.07, 0.13, 0.15],
                            'H': [2.99, 0.23, 4.66, 0.13, 7.69, 0.27, 0.30],
                            'I': [4.19, 0.19, 4.00, 1.80, 6.04, 0.30, 0.45],
                            'K': [1.89, 0.22, 4.77, -0.99, 9.99, 0.32, 0.27],
                            'L': [2.59, 0.19, 4.00, 1.70, 6.04, 0.39, 0.31],
                            'M': [2.35, 0.22, 4.43, 1.23, 5.71, 0.38, 0.32],
                            'N': [1.60, 0.13, 2.95, -0.60, 6.52, 0.21, 0.22],
                            'P': [2.67, 0.00, 2.72, 0.72, 6.80, 0.13, 0.34],
                            'Q': [1.56, 0.18, 3.95, -0.22, 5.65, 0.36, 0.25],
                            'R': [2.34, 0.29, 6.13, -1.01, 10.74, 0.36, 0.25],
                            'S': [1.31, 0.06, 1.60, -0.04, 5.70, 0.20, 0.28],
                            'T': [3.03, 0.11, 2.60, 0.26, 5.60, 0.21, 0.36],
                            'V': [3.67, 0.14, 3.00, 1.22, 6.02, 0.27, 0.49],
                            'W': [3.21, 0.41, 8.08, 2.25, 5.94, 0.32, 0.42],
                            'Y': [2.94, 0.30, 6.47, 0.96, 5.66, 0.25, 0.41]}

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    pdb_name=sequence_name
    #print pdb_name
    sequence_file=open(input_dir,'r').readlines()
    output_file=open(output_dir+'/'+pdb_name+'.7PC','w')

    phys = [phys_dic.get(i, phys_dic['A']) for i in sequence_file[1].strip()]

    input_feature = window_data(phys)

    for i in xrange(len(input_feature)):
        for j in xrange(len(input_feature[i])):
             output_file.write(str(input_feature[i][j]) + ' ')
        output_file.write('\n')
    output_file.close()

    pdb_name=sequence_name
    #print pdb_name
    normalize_7PC(output_dir+'/'+pdb_name+'.7PC',final_output_dir,pdb_name,max_7PC,min_7PC)

#=====================================finish 7PC=======================================================================

#=====================================get disorder=====================================================================
def get_input_DISORDER(disorder_file):
    fp = open(disorder_file, 'r')
    lines = fp.readlines()
    fp.close()

    disorder = []
    
    for i in xrange(len(lines[2].strip().split(' '))):
        tmp_disorder=[]
        tmp_disorder.append(float(lines[2].strip().split(' ')[i]))
        disorder.append(tmp_disorder)

    return disorder


def normalize_DISORDER(input_disorder,final_output,min_max_dir,pdb_name,max_disorder,min_disorder):
    input_dir=input_disorder
    output_dir=final_output

    max_list=open(max_disorder,'r').readlines()
    min_list=open(min_disorder,'r').readlines()

    feature_list=open(input_dir,'r').readlines()
    output_feature_file=open(output_dir+'/'+pdb_name+'_disorder','w')
    for i in xrange(len(feature_list)):

        for j in xrange(len(feature_list[i].split(' '))-1):
            feature_value=float(float(feature_list[i].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))/(float(max_list[0].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))
            if(float(feature_value)>1):
                feature_value=1.0
            if(float(feature_value)<0.0):
                feature_value=0.0
            output_feature_file.write(str(feature_value)+' ')
        output_feature_file.write('\n')
    output_feature_file.close()

def get_DISORDER(feature_number,normalize,in_dir,out_dir,final_output_dir,winsize,normalize_dir,sequence_name):

    if_normalize=normalize
    input_dir=in_dir
    output_dir=out_dir

    max_disorder=normalize_dir+'/max_disorder'
    min_disorder=normalize_dir+'/min_disorder'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    pdb_name=sequence_name
    #print pdb_name
    disorder_file=input_dir
    output_file=open(output_dir+'/'+pdb_name+'.DISORDER','w')

    disorder = get_input_DISORDER(disorder_file)
    #print disorder
    input_feature = window_data(disorder)

    for i in xrange(len(input_feature)):
        for j in xrange(len(input_feature[i])):
             output_file.write(str(input_feature[i][j]) + ' ')
        output_file.write('\n')
    output_file.close()

    #normalize DISORDER
    pdb_name=sequence_name
    #print pdb_name
    normalize_DISORDER(output_dir+'/'+pdb_name+'.DISORDER',final_output_dir,out_dir,pdb_name,max_disorder,min_disorder)

#==================================finish disorder=====================================================================

#============================get contacts=============================================================================

def get_input_CONTACTS_15(contacts_file,feature_number):
    #open the two files, read in their data and then close them
    fp = open(contacts_file, 'r')
    lines = fp.readlines()
    fp.close()

    contacts=[]

    for line in lines:
        split_line=line.split()
        contacts_tmp=[]
        if(len(split_line) == int(feature_number)):
            for j in xrange(len(split_line)):
                contacts_tmp.append(split_line[j])    
        contacts.append(contacts_tmp)
    return contacts

def normalize_CONTACTS_15(input_contacts,final_output,min_max_dir,pdb_name,max_contacts15,min_contacts15):
    input_dir=input_contacts

    output_dir=final_output
    input_max=max_contacts15
    input_min=min_contacts15

    max_list=open(input_max,'r').readlines()
    min_list=open(input_min,'r').readlines()

    feature_list=open(input_dir,'r').readlines()
    output_feature_file=open(output_dir+'/'+pdb_name+'_contacts_15','w')
    for i in xrange(len(feature_list)):
        for j in xrange(len(feature_list[i].split(' '))-1):
            feature_value=float(float(feature_list[i].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))/(float(max_list[0].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))
            if(float(feature_value)>1):
                feature_value=1.0
            if(float(feature_value)<0.0):
                feature_value=0.0
            output_feature_file.write(str(feature_value)+' ')
        output_feature_file.write('\n')

    output_feature_file.close()

def get_CONTACTS_15(feature_number,normalize,in_dir,out_dir,final_output_dir,winsize,normalize_dir,sequence_name):

    if_normalize=normalize
    input_dir=in_dir
    output_dir=out_dir
    max_contacts15=normalize_dir+'/max_contacts15'
    min_contacts15=normalize_dir+'/min_contacts15'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    pdb_name=sequence_name
    #print pdb_name
    contacts_file=input_dir
    output_file=open(output_dir+'/'+pdb_name+'.CONTACTS_15','w')

    contacts = get_input_CONTACTS_15(contacts_file,feature_number)
        
    input_feature = window_data(contacts)

    for i in xrange(len(input_feature)):
        for j in xrange(len(input_feature[i])):
             output_file.write(str(input_feature[i][j]) + ' ')
        output_file.write('\n')
    output_file.close()

    #normalize CONTACTS
    db_name=sequence_name
    #print pdb_name
    normalize_CONTACTS_15(output_dir+'/'+pdb_name+'.CONTACTS_15',final_output_dir,out_dir,pdb_name,max_contacts15,min_contacts15)

#============================ finish contacts==========================================================================

#=============================get fragments============================================================================
def get_input_fragments(fragments_file,feature_number,mean_Dev_list,pdb_name,output_dir,sequence_file):

    fp = open(fragments_file, 'r')
    lines = fp.readlines()
    fp.close()

    # declare the empty dictionary with each of the entries
    #sequence_dir='/exports/store2/haiou/JieShare/Project_of_TorsionAngle/Inputs'
    sequence_file=open(sequence_file,'r').readlines()
    fragments=[]

    new_frag_out = open(output_dir+'/'+pdb_name+'_new_frag', 'w')
    tmp_i=0
    for line in lines:
        split_line=line.split()

        if(len(split_line) == (int(feature_number) + 1)  and (split_line[0] != 'ID')):
            fragments_tmp= [float(i) for i in split_line[0:int(feature_number)]]
            for j in xrange(len(mean_Dev_list)):
                if(mean_Dev_list[j][0] == sequence_file[1][tmp_i]):
                    mean_phi=mean_Dev_list[j][1]
                    var_phi=mean_Dev_list[j][2]
                    mean_psi=mean_Dev_list[j][3]
                    var_psi=mean_Dev_list[j][4]
            tmp_i+=1

            tmp_phi_middle=float(fragments_tmp[1])-float(mean_phi)
            tmp_phi_std=float(var_phi)

            tmp_psi_middle=float(fragments_tmp[3])-float(mean_psi)
            tmp_psi_std=float(var_psi)

            if(float(tmp_phi_middle)>180):
                tmp_phi_middle=180
            elif(float(tmp_phi_middle)<-180):
                tmp_phi_middle=-180

            if(float(tmp_phi_std)>180):
                tmp_phi_std=180
            elif(float(tmp_phi_std)<-180):
                tmp_phi_std=-180

            if(float(tmp_psi_middle)>180):
                tmp_psi_middle=180
            elif(float(tmp_psi_middle)<-180):
                tmp_psi_middle=-180

            if(float(tmp_psi_std)>180):
                tmp_psi_std=180
            elif(float(tmp_psi_std)<-180):
                tmp_psi_std=-180

            phi_middle=math.radians(tmp_phi_middle)
            phi_right=math.radians(tmp_phi_std)

            psi_middle=math.radians(tmp_psi_middle)
            psi_right=math.radians(tmp_psi_std)
            
            tmp_=[]
            tmp_.append(phi_middle)
            tmp_.append(phi_right)

            tmp_.append(psi_middle)
            tmp_.append(psi_right)

            fragments.append(tmp_)

    new_frag_out.close()

    return fragments

def normalize_FRAGMENTS_error(input_fragments,final_output,min_max_dir,pdb_name,max_fragments,min_fragments):
    input_dir=input_fragments
    output_dir=final_output

    max_list=open(max_fragments,'r').readlines()
    min_list=open(min_fragments,'r').readlines()

    feature_list=open(input_dir,'r').readlines()
    output_feature_file=open(output_dir+'/'+pdb_name+'_fragments_error','w')
    for i in xrange(len(feature_list)):
        for j in xrange(len(feature_list[i].split(' '))-1):
            feature_value=float(float(feature_list[i].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))/(float(max_list[0].split(' ')[j].strip())-float(min_list[0].split(' ')[j].strip()))
            if(float(feature_value)>1):
                feature_value=1.0
            if(float(feature_value)<0.0):
                feature_value=0.0
            output_feature_file.write(str(feature_value)+' ')
        output_feature_file.write('\n')
    output_feature_file.close()

def get_FRAGMENTS_error(feature_number,normalize,in_dir,out_dir,final_output_dir,winsize,normalize_dir,sequence_name,sequence_file):

    if_normalize=normalize
    input_dir=in_dir
    output_dir=out_dir

    max_fragments=normalize_dir+'/max_fragments'
    min_fragments=normalize_dir+'/min_fragments'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    mean_Dev_list=[
    ['A', '2.33813826398852', '32.6375986891078', '3.87284667503587', '61.7313401916215'],
    ['C', '0', '0', '0', '0'],
    ['D', '0', '0', '0', '0'],
    ['E', '4.71041226388889', '34.1964773619604', '1.35061577314815', '60.8252265527372'],
    ['F', '8.60854856240126', '32.386575912434', '-0.465576642969983', '59.5907448059746'],
    ['G', '0', '0', '0', '0'],
    ['H', '0', '0', '0', '0'],
    ['I', '3.38362404419322', '28.0036014066794', '5.72439400308325', '56.8805622191122'],
    ['K', '3.81206886726547', '40.9159358704877', '0.988647135728543', '63.6288064078162'],
    ['L', '4.23341747855392', '29.2642997732282', '-2.37553730392157', '57.7526756531369'],
    ['M', '-0.823469890965731', '40.2352546300108', '-2.06854317757009', '65.3715772972374'],
    ['N', '0', '0', '0', '0'],
    ['P', '-8.0343925990099', '19.8133022240976', '-11.1551639542079', '81.3386859470953'],
    ['Q', '5.28255321376281', '38.0316654002613', '2.73349969253294', '62.0922662452154'],
    ['R', '4.38199269430052', '39.2161784211576', '0.175437761658032', '62.8262138932591'],
    ['S', '0', '0', '0', '0'],
    ['T', '0', '0', '0', '0'],
    ['V', '2.27951570209465', '28.8798641053739', '4.71833349883631', '57.5389159493841'],
    ['W', '4.50126837254902', '34.8395628025976', '-4.69090215686274', '59.1260827231632'],
    ['Y', '8.5099081795302',  '34.692968591627', '-2.63551320469799', '60.7086535874481']
    ]
    
    #======================read dataset on whole dataset===================================================
    pdb_name=sequence_name
    #print pdb_name
    fragments_file=input_dir
    output_file=open(output_dir+'/'+pdb_name+'.tmp_FRAGMENTS_error','w')

    fragments = get_input_fragments(fragments_file,feature_number,mean_Dev_list,pdb_name,output_dir,sequence_file)
    input_feature = window_data(fragments)

    for i in xrange(len(input_feature)):
        for j in xrange(len(input_feature[i])):
             output_file.write(str(input_feature[i][j]) + ' ')
        output_file.write('\n')
    output_file.close()

    #normalize FRAGMENTS on training dataset
    pdb_name=sequence_name
    normalize_FRAGMENTS_error(output_dir+'/'+pdb_name+'.tmp_FRAGMENTS_error',final_output_dir,out_dir,pdb_name,max_fragments,min_fragments)    
    #print 'finish fragments error normalize'

#================================finish fragments========================================================================

#==================================combine features======================================================================
def combine_features_and_number(feature_list,winsize,final_output_dir,sequence_name):
    window=int(winsize)

    name=''
    for i in xrange(len(feature_list)):
        name=name+str(feature_list[i][0])+str('_')

    feature_combine_output=final_output_dir

    pdb_name=sequence_name
    #print pdb_name

    whole_list_for_each_sequence=[]

    total_column=0
    for i in xrange(len(feature_list)):
        feature_number=feature_list[i][1]

        total_column=total_column+int(feature_number)

        if(feature_list[i][0] == 'our_new_PSSM' and feature_list[i][1] == '20'):
            input_pssm_file=open(final_output_dir+'/'+pdb_name+'_pssm_own_new','r').readlines()
            whole_list_for_each_sequence.append(input_pssm_file)
            #print "add our_new_PSSM\n"

        if(feature_list[i][0] == '7PC' and feature_list[i][1] == '7'):
            input_7pp_file=open(final_output_dir+'/'+pdb_name+'_7pc','r').readlines()
            whole_list_for_each_sequence.append(input_7pp_file)
            #print 'add 7PP\n'

        if(feature_list[i][0] == 'SS8' and feature_list[i][1] == '8'):
            input_8ss_file=open(final_output_dir+'/'+pdb_name+'_8_ss','r').readlines()
            whole_list_for_each_sequence.append(input_8ss_file)
            #print 'add SS8\n'

        if(feature_list[i][0] == 'ASA' and feature_list[i][1] == '1'):
            input_1asa_file=open(final_output_dir+'/'+pdb_name+'_1_asa','r').readlines()
            whole_list_for_each_sequence.append(input_1asa_file)
            #print 'add ASA\n'

        if(feature_list[i][0] == 'DISORDER' and feature_list[i][1] == '1'):
            input_disorder_file=open(final_output_dir+'/'+pdb_name+'_disorder','r').readlines()
            whole_list_for_each_sequence.append(input_disorder_file)
            #print 'add DISORDER\n'

        if(feature_list[i][0] == 'CONTACTS15' and feature_list[i][1] == '15'):
            input_disorder_file=open(final_output_dir+'/'+pdb_name+'_contacts_15','r').readlines()
            whole_list_for_each_sequence.append(input_disorder_file)
            #print "add CONTACTS15\n"

        if(feature_list[i][0] == 'FRAGMENTS_error' and feature_list[i][1] == '4'):
            input_disorder_file=open(final_output_dir+'/'+pdb_name+'_fragments_error','r').readlines()
            whole_list_for_each_sequence.append(input_disorder_file)

        #============================================renumber the features================
        #print "total_column===="+str(total_column)
        total_column_each_R=int(total_column)*int(window)
        output_file=open(feature_combine_output+'/'+pdb_name+'.fea_'+str(winsize),'w')
        for i in xrange(len(whole_list_for_each_sequence[0])): 
            flag=1

            for j in xrange(len(whole_list_for_each_sequence)):#5
                for k in xrange(len(whole_list_for_each_sequence[j][i].split())):
                    output_file.write(str(flag)+":"+str(whole_list_for_each_sequence[j][i].split()[k])+' ')
                    flag=flag+1
            output_file.write('\n')
        output_file.close()
    
#=========================================================================================================================

if __name__ == '__main__':

    if len(sys.argv) != 6:
        print '\n'
        print 'please input the right parametets,sequence_name, outputdir, install_dir, sequence_file, modeltype '
        sys.exit(1)

    sequence_name=sys.argv[1]
    final_output_dir=sys.argv[2]
    install_dir=sys.argv[3]
    sequence_file=sys.argv[4]
    modeltype=sys.argv[5]
    
    
    # set modeltype and window size
    set_windowsize_by_modeltype(modeltype)

    normalize_dir=sys.argv[3]+'/Inputs/normalize/win'+str(winsize)

    start_time = time.time()
    
    feature_list=[

                    ['PSSM','20','1',final_output_dir+'/features/tmp/'+sequence_name+'.pssm',
                    final_output_dir+'/features/tmp/'],    

                    ['ASA','1','0',final_output_dir+'/features/tmp/'+sequence_name+'.ss_ASA.acc',
                    final_output_dir+'/features/tmp/'],

                    ['SS8','8','0',final_output_dir+'/features/tmp/'+sequence_name+'.ss_ASA.ss8',
                    final_output_dir+'/features/tmp/'],

                    ['7PC','7','1',sequence_file,
                    final_output_dir+'/features/tmp/'],

                    ['DISORDER','1','1',final_output_dir+'/features/tmp/'+sequence_name+'.disorder',
                    final_output_dir+'/features/tmp/'],
        
                    ['CONTACTS_15','15','1',final_output_dir+'/features/tmp/'+sequence_name+'.contacts15',
                    final_output_dir+'/features/tmp/'],    
                
                    ['FRAGMENTS_error','4','1',final_output_dir+'/features/tmp/'+sequence_name+'.fragments',
                    final_output_dir+'/features/tmp/'],                
    ]    
    
    
    for i in xrange(len(feature_list)):
        if(feature_list[i][0] == 'PSSM'):
            get_PSSM(feature_list[i][1],feature_list[i][2],feature_list[i][3],feature_list[i][4],normalize_dir,final_output_dir,winsize,sequence_name)
        
        if(feature_list[i][0] == 'ASA'):
            get_ASA(feature_list[i][1],feature_list[i][2],feature_list[i][3],feature_list[i][4],final_output_dir,winsize,sequence_name)

        if(feature_list[i][0] == 'SS8' and feature_list[i][1] == '8'):
            get_SS_8class(feature_list[i][1],feature_list[i][2],feature_list[i][3],feature_list[i][4],final_output_dir,winsize,sequence_name)

        if(feature_list[i][0] == '7PC'):
            get_7PC(feature_list[i][1],feature_list[i][2],feature_list[i][3],feature_list[i][4],final_output_dir,winsize,normalize_dir,sequence_name)
    
        if(feature_list[i][0] == 'DISORDER'):
            get_DISORDER(feature_list[i][1],feature_list[i][2],feature_list[i][3],feature_list[i][4],final_output_dir,winsize,normalize_dir,sequence_name)

        if(feature_list[i][0] == 'CONTACTS_15'):
            get_CONTACTS_15(feature_list[i][1],feature_list[i][2],feature_list[i][3],feature_list[i][4],final_output_dir,winsize,normalize_dir,sequence_name)

        if(feature_list[i][0] == 'FRAGMENTS_error'):
            get_FRAGMENTS_error(feature_list[i][1],feature_list[i][2],feature_list[i][3],feature_list[i][4],final_output_dir,winsize,normalize_dir,sequence_name,sequence_file)


    #['our_new_PSSM','20'],['SS8','8'],['7PC','7'],['CONTACTS15','15'],['DISORDER','1'],['ASA','1'],['FRAGMENTS_error','4']
    #['our_new_PSSM','20'],['SS8','8'],['7PC','7'],['CONTACTS15','15'],['DISORDER','1'],['ASA','1'],['FRAGMENTS_error','4']
    #old ['our_new_PSSM','20'],['SS8','8'],['CONTACTS15','15'],['DISORDER','1'],['FRAGMENTS_error','4'],['7PC','7'],['ASA','1']
    combine_features_and_number([['our_new_PSSM','20'],['SS8','8'],['7PC','7'],['CONTACTS15','15'],['DISORDER','1'],['ASA','1'],['FRAGMENTS_error','4']],
        winsize,final_output_dir,sequence_name)

    end_time = time.time()
    
    print "# The features generation time is: "+str(float(end_time-start_time)/60.0)+'min'
    print "# The features file is saved to:   "+final_output_dir+'/'+sequence_name+'.fea_'+str(winsize)
    print "###############################################################################\n"
