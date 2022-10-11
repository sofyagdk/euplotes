
import os
for i in range(32):
	nowfile = 'gra_' + str(i) + '.sh'
	with open(nowfile, 'w') as fout:
		fout.write("#!/bin/bash\n#$ -cwd\npython graphmaker_onlylog_v5.py -f '")
		fout.write(str((i+1) * 50))
		fout.write("' -s " + str((i) * 50) + " -e ")
		fout.write(str((i+1) * 50))
	os.system("qsub gra_" + str(i) + '.sh')


nowfile = 'gra_33.sh'
with open(nowfile, 'w') as fout:
	fout.write("#!/bin/bash\n#$ -cwd\npython graphmaker_onlylog_v5.py -f '")
	fout.write(str(1634))
	fout.write("' -s " + str(1600) + " -e ")
	fout.write(str(1634))
os.system("qsub gra_33.sh")
