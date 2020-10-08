#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.beams3d import read_beams3d

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

lfail = 0
failtol = 5.0
filename='beams3d_ORBITS_dist.h5'
data=read_beams3d(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

# Calc values

print('BEAMS3D VERSION: ' + str(round(data['VERSION'],2)))
print('==== Matrix ====')
varlist={}
varlist['dense_prof']=np.array([6.76E12, 4.32E12, 5.56E12, 5.42E12])
varlist['ipower_prof']=np.array([1.42, 0.90, 0.96, 0.93])
varlist['epower_prof']=np.array([2.47, 1.58, 2.35, 2.30])
varlist['j_prof']=np.array([-0.022, -0.068, 0.0057, -0.092])
#print(data.keys())
#print(data['dense_prof'][16,0])
#print(data['dense_prof'][16,1])
#print(data['dense_prof'][32,0])
#print(data['dense_prof'][32,1])
for temp in varlist:
    ir=16; ib = 0;
    act = varlist[temp][0]
    cal = data[temp][ir,ib]
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    #perct = max(abs(act-cal))
    #print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(int(perct))+'%')
    perct = 100*abs((act-cal)/act)
    print('  '+temp+'['+str(ir)+','+str(ib)+']: '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
    ir=16; ib = 1;
    act = varlist[temp][1]
    cal = data[temp][ir,ib]
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    #perct = max(abs(act-cal))
    #print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(int(perct))+'%')
    perct = 100*abs((act-cal)/act)
    print('  '+temp+'['+str(ir)+','+str(ib)+']: '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
    ir=32; ib = 0;
    act = varlist[temp][2]
    cal = data[temp][ir,ib]
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    #perct = max(abs(act-cal))
    #print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(int(perct))+'%')
    perct = 100*abs((act-cal)/act)
    print('  '+temp+'['+str(ir)+','+str(ib)+']: '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
    ir=32; ib = 1;
    act = varlist[temp][3]
    cal = data[temp][ir,ib]
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    #perct = max(abs(act-cal))
    #print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(int(perct))+'%')
    perct = 100*abs((act-cal)/act)
    print('  '+temp+'['+str(ir)+','+str(ib)+']: '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




