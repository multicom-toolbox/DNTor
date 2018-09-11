import sys
import os
import os.path
import math
import shutil
import numpy 
import time
import subprocess
import errno  


def generate_PSSM(seq_name,seq_file,output_dir,pssm_tool_dir):

	exact_script=pssm_tool_dir+'/script/generate_flatblast.pl'
	script_dir=pssm_tool_dir+'/script/'
	blast_dir=pssm_tool_dir+'/blast-2.2.17/bin'
	big_db=pssm_tool_dir+'/db/big/big_98_X'
	nr_db=nr_database_dir+'/nr'
	#seq_file=
	output_prefix_alg=output_dir+seq_name


	#print 'perl script_dir'+' blast_dir script_dir big_db nr_db seq_file output_prefix_alg >/home/lihaio/JieShare/ProteinTorsionAngle_Prediction/Outputs/features/tmp/tmp '
	subprocess.call(["perl", exact_script, blast_dir, script_dir, big_db, nr_db, seq_file, output_prefix_alg ])
	#+'>/home/lihaio/JieShare/ProteinTorsionAngle_Prediction/Outputs/features/tmp/tmp'


def generate_Disorder(seq_name,seq_file,output_dir,disorder_tool_dir):
	exact_script=disorder_tool_dir+'/bin/predict_diso.sh'

	os.system(exact_script+' '+seq_file+' '+output_dir+seq_name+'.disorder')

def generate_SS_ASA(seq_name,seq_file,output_dir,SCRATCH_tool_dir):
	exact_script=SCRATCH_tool_dir+'/bin/run_SCRATCH-1D_predictors.sh'
	os.system(exact_script+' '+seq_file+' '+output_dir+seq_name+'.ss_ASA 4')

def generate_Contacts(seq_name,seq_file,output_dir,contacts_tool_dir):
	exact_script_for_gtg=contacts_tool_dir+'buildFeature'
	exact_script_for_Pred=contacts_tool_dir+'AcconPred'	
	#change the directory
	os.chdir(contacts_tool_dir)
	os.system('pwd')
	#print exact_script_for_gtg+' -i '+seq_file+' -o '+output_dir+seq_name+'.tgt -c 4'
	os.system(exact_script_for_gtg+' -i '+seq_file+' -o '+output_dir+seq_name+'.tgt -c 4')
	os.chdir(output_dir)
	#os.system('pwd')
	os.system(exact_script_for_Pred+' '+output_dir+seq_name+'.tgt 0 >'+output_dir+seq_name+'.tmp_contacts15')

	input_contacts_file=open(output_dir+seq_name+'.tmp_contacts15','r').readlines()

	output_contacts_file=open(output_dir+seq_name+'.contacts15','w')

	for i in xrange(5,len(input_contacts_file)):
		for j in xrange(4,19):
			tmp_contact=input_contacts_file[i].split(' ')[j]
			output_contacts_file.write(str(tmp_contact)+' ')
		output_contacts_file.write('\n')
	output_contacts_file.close()
	os.system('rm '+contacts_tool_dir+ '/tmp/'+seq_name+'*')

def generate_Sspro(seq_name,seq_file,output_dir,sspro_tool_dir):
	exact_script=sspro_tool_dir+'/bin/predict_ssa_ab.sh'
	os.system(exact_script+' '+seq_file+' '+output_dir+seq_name+'.sspro')
	origin_ss=open(output_dir+seq_name+'.sspro','r').readlines()
	final_output=open(output_dir+seq_name+'.ss4fragments','w')
	final_output.write('>')
	final_output.write(origin_ss[0])
	final_output.write(origin_ss[2])

	final_output.close()

def generate_Fragments(seq_name,seq_file,output_dir,fragments_tool_dir):
        exact_script=fragments_tool_dir+'/FRAGSION_linux64'
        model=fragments_tool_dir+'/IOHMM.dat'
        ss=output_dir+'/'+seq_name+'.ss4fragments'
        output_file=open(output_dir+seq_name+'.fragments','w')
        output_file.write('ID   PHI_mean        PHI_std PSI_mean        PSI_std\n')

        os.system(exact_script+' -f '+seq_file+' -s '+ss+' -m '+model+' -o '+output_dir+seq_name+'.tmp_fragments')

        #extract fragments
        sequence_len=len(open(seq_file,'r').readlines()[1].strip())
        fra_file=open(output_dir+seq_name+'.tmp_fragments', 'r').readlines()

        one_AA_frag_list=[]
        #print len(fra_file)
        for i in xrange(2,len(fra_file),802):
                #print i
                if(fra_file[i].split()[0] == 'xxxx'):
                        tmp_frag_list=[]
                        i=i+2
                        for j in xrange(0,200):
                                #print i
                                tmp_frag_list.append(fra_file[i].split()[5] + ' '+fra_file[i].split()[6])
                                i=i+4
                                #print tmp_frag_list
                        one_AA_frag_list.append(tmp_frag_list)

        tmp_frag_list_last2=[]
        #for i in xrange(68173,len(fra_file),4):
        for i in xrange(len(fra_file)-800+2,len(fra_file),4):
        	#print len(fra_file)-800+2
                tmp_frag_list_last2.append(fra_file[i].split()[5] + ' '+fra_file[i].split()[6])
        one_AA_frag_list.append(tmp_frag_list_last2)

        tmp_frag_list_last1=[]
        #for i in xrange(68174,len(fra_file),4):
        for i in xrange(len(fra_file)-800+3,len(fra_file),4):
                tmp_frag_list_last1.append(fra_file[i].split()[5] + ' '+fra_file[i].split()[6])
        one_AA_frag_list.append(tmp_frag_list_last1)

        #print len(one_AA_frag_list)

        for i in xrange(len(one_AA_frag_list)):
                #print one_AA_frag_list[i]
                phi_list=[]
                psi_list=[]

                for j in xrange(len(one_AA_frag_list[i])):
                        #print one_AA_frag_list[i][j].split()[0]
                        phi_list.append(float(one_AA_frag_list[i][j].split()[0]))
                        psi_list.append(float(one_AA_frag_list[i][j].split()[1]))

                phi_mean=numpy.mean(phi_list)
                psi_mean=numpy.mean(psi_list)

                phi_std=numpy.std(phi_list)
                psi_std=numpy.std(psi_list)

                output_file.write(str(i+1)+'    '+str(phi_mean)+'       '+str(phi_std)+'        '+str(psi_mean)+'       '+str(psi_std)+'\n')

        output_file.close()

        #os.system('rm '+fragments_tool_dir+ '/tmp/'+seq_name+'*')
        #FRAGSION -f fasta -s ss -m model -o outfile


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

if __name__ == '__main__':

	#print len(sys.argv)
	if len(sys.argv) != 5:
		print 'please input the right parameters sequence_name,sequence_file,install_dir, output_dir'
		sys.exit(1)

	sequence=sys.argv[1]
	sequence_file=sys.argv[2]
	install_dir=sys.argv[3]
	outdir=sys.argv[4]

	input_file=sequence_file
	output_dir=outdir+'/features/tmp/'
	if os.path.isdir(output_dir):
		print "Empty feature folder: ",output_dir
		filelist = [ f for f in os.listdir(output_dir)]
		for f in filelist:
			print "remove file: ",f
			file_to_remove = output_dir +'/'+f
			os.remove(file_to_remove)
	else:
		print "Creating feature folder: ",output_dir
		mkdir_p(output_dir)
	
	nr_database_dir=install_dir+'/Methods/Tools/database/nr/'
	
	#produce PSSM
	pssm_tool_dir=install_dir+'/Methods/Tools/pspro2.0/' 
	generate_PSSM(seq_name=sequence,seq_file=input_file,output_dir=output_dir,pssm_tool_dir=pssm_tool_dir)

	#produce Disorder
	disorder_tool_dir=install_dir+'/Methods/Tools/predisorder1.1/'
	generate_Disorder(seq_name=sequence,seq_file=input_file,output_dir=output_dir,disorder_tool_dir=disorder_tool_dir)

	#SS and ASA
	SCRATCH_tool_dir=install_dir+'/Methods/Tools/SCRATCH-1D_1.1/'
	generate_SS_ASA(seq_name=sequence,seq_file=input_file,output_dir=output_dir,SCRATCH_tool_dir=SCRATCH_tool_dir)

	#Contacts
	contacts_tool_dir=install_dir+'/Methods/Tools/AcconPred_package_v1.00/'
	generate_Contacts(seq_name=sequence,seq_file=input_file,output_dir=output_dir,contacts_tool_dir=contacts_tool_dir)

	#SS for fragments
	print "sspro running"
	sspro_tool_dir=install_dir+'/Methods/Tools/sspro4/'
	generate_Sspro(seq_name=sequence,seq_file=input_file,output_dir=output_dir,sspro_tool_dir=sspro_tool_dir)

	#Fragments
	fragments_tool_dir=install_dir+'/Methods/Tools/fragsion/'
	generate_Fragments(seq_name=sequence,seq_file=input_file,output_dir=output_dir,fragments_tool_dir=fragments_tool_dir)

