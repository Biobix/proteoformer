import traceback
import os
import sys
import re
import mpl_toolkits.mplot3d
location = mpl_toolkits.mplot3d.__path__[0]
pyfile = location+"/axes3d.py"
backup = location+"/axes3d_old.py"
pycfile = location+"/axes3d.pyc"

print "\nLocation of the mplot3d file that should have been modified: "+pyfile+"\n"

found = False
input = open(pyfile, 'r')
for line in input:
	m = re.search(re.escape("# Small modifications for cuboidal 3D plotting area <steven.verbruggen@ugent.be>"), line)
	if m:
		found = True
		break
input.close()

if found==True:
	print "File already modified\n"
else:
	os.system("mv "+pyfile+" "+backup)
	os.system("rm -rf "+pycfile)
	os.system("wget -q --no-check-certificate https://raw.githubusercontent.com/Biobix/mQC/master/mqc_tools/site-packages/axes3d.py")
	os.system("mv axes3d.py "+location)
	print "Axes3d file modified\n"
