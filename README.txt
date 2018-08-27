#########################################################################################################################
DeepQA: Master estimation of model accuracy with deep neural networks
Free for Academic only.
All rights reserved.

Authors:
Renzhi Cao             rcrg4@mail.missouri.edu
debswapna bhattacharya db279@mail.missouri.edu
Jie Hou                jh7x3@mail.missouri.edu
Jianlin Cheng          chengji@missouri.edu         
Department of Computer Science
University of Missouri, Columbia

Please send your questions, feedback and comments to: chengji@missouri.edu.

#########################################################################################################################
1. INSTALLATION
#########################################################################################################################

The software was developed and tested on the RedHat 4.4.7 and ARCH Linux 3.8.11

1). Download the DeepQA tool (DeepQA.tar.gz) at http://sysbio.rnet.missouri.edu/multicom_toolbox/tools.html#license, 
	The DeepQA tool includes:
	--------------------------------------------
	Name                Size          Size
					(decompressed) (compressed)
	--------------------------------------------
	README.txt      4.00   K      
	bin/            8.00   K      
	data/           4.00   K      
	script/         392    K      
	test/           132    K      

	TOTAL          ~5.80   GB        ~2.00  GB
	--------------------------------------------
2). Unzip DeepQA_pacage.tar.gz
	$ tar -zxvf DeepQA_package.tar.gz

3). Go to the DeepQA folder, run the configure file
	$ ./configure.pl

4). Installation is done if you don't see any errors, otherwise fix it by installing the missing tools or contact the authors.

-------------------------------------------------------------------------------------------------------------------------

There should be one file named DeepQA.sh in the bin directory, simply run it and you will see three inputs are needed. path of fasta sequence, directory of input pdbs, directory of output.

(1). Path of sequence in fasta format. The file is in fasta format, and the sequence inside should be the sequence of the 3D model (pdb) to be evaluated.

(2). Path of the folder including all models for evaluation. All models in pdb format with the same sequence is put together inside this input folder. 

(3). Path of the output folder. The evaluation result (DeepQA.scores) will be stored in this folder. 

     For the final score, the first column in the output prediction file is the model name, the second column is the global quality score. 


-------------------------------------------------------------------------------------------------------------------------

#########################################################################################################################
2. USAGE
#########################################################################################################################

Testing example:
1) Go to test folder in the DeepQA tool, run:
	$ cd ./test/
	$ ../bin/DeepQA.sh T0709.fasta T0709 test_T0709

   The predicted score can be found in 'test_T0709/Predictions.txt'

-------------------------------------------------------------------------------------------------------------------------
