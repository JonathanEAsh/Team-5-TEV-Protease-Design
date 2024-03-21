#from pyrosetta import *
#init()
import os
import os.path
import time
#from joey_utils import index_selector
#from joey_utils import intergroup_selector
#from joey_utils import make_task_factory
#from joey_utils import fast_relax_mover
#from joey_utils import extract_pose_chain
import subprocess
import pandas as pd
# df=pd.read_csv("incompletes.csv")
# col=df['ID']
# names=[]
# for x in col:
# 	names.append(x)
#with open('test.fa') as f:
with open('final_seqs_meth.fa') as f:
	ctr=1
	design_ctr=0
	for line in f:
		if ctr%2==0:
			design_ctr+=1
			if os.path.isfile('relax_output/'+temp_name+'.txt')==False:
			#if full_name in names:
				print(temp_name)
				design_seq=line.strip()
				process1=subprocess.Popen(['squeue','-u','ja961'], stdout=subprocess.PIPE)
				process2=subprocess.check_output(['wc','-l'], stdin=process1.stdout)
				while int(process2)>=500:
					time.sleep(300)
					process1=subprocess.Popen(['squeue','-u','ja961'], stdout=subprocess.PIPE)
					process2=subprocess.check_output(['wc','-l'], stdin=process1.stdout)
				#command="sbatch python relax_run.py --og_seq="+og_seq+" --design_seq="+design_seq+" --design_num="+str(design_ctr)
				subprocess.run(['sbatch', '--time=2:00:00','--mem=32G',"./relax_submit.sh",str(temp_name), design_seq])
				#'--error=/dev/null','./relax_submit.sh', 
				#subprocess.run(['python', 'slurmit_BAY.py', '--job=relax_'+str(design_ctr) --partition main --tasks 1 --cpus 1 --mem 16G --time 20:00:00 --begin now --requeue True --outfiles $od/PDZ_structures/job_${jobnum}/logs/${jobname}_%a --command "$com"])
				#subprocess.run(["sbatch", "python", "relax_run.py", "--og_seq="+og_seq, "--design_seq="+design_seq, "--design_num="+str(design_ctr)])
		else:
			temp_name=line.replace('>','').strip()
		ctr+=1
#python slurmit_BAY.py --job $jobname --partition main --tasks 1 --cpus 1 --mem 16G --time 20:00:00 --begin now --requeue True --outfiles $od/PDZ_structures/job_${jobnum}/logs/${jobname}_%a --command "$com"
