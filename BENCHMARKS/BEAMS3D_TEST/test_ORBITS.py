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
failtol = 1.0
filename='beams3d_ORBITS.h5'
data=read_beams3d(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)

# Calc values
rho = np.sqrt(data['S_lines'])
rho_max = np.max(rho,axis=1)
rho_min = np.min(rho,axis=1)
data['delta'] = rho_max-rho_min
x = data['R_lines']-10.0
y = data['Z_lines']
theta = np.arctan2(y,x)
theta = np.where(theta > np.pi,theta-pi,theta)
data['turning'] = np.max(theta,axis=1)
data['R0']=data['R_lines'][:,0]
data['R1']=data['R_lines'][:,1]
data['R100']=data['R_lines'][:,100]
data['R500']=data['R_lines'][:,500]

print('BEAMS3D VERSION: ' + str(data['VERSION']))
print('==== Vectors ====')
varlist={}
varlist['turning']=np.array([0.12275448, 0.1788069,  0.23483922, 0.29094046, 0.34722868, 0.40375728,
 0.46065893, 0.51799329, 0.57587523, 0.63441059, 0.69368621, 0.75391989,
 0.81514834, 0.87755329, 0.94131251, 1.0066042,  1.073641,   1.14261056,
 1.21404663, 1.2880335,  1.36508321, 1.44570684, 1.53053939, 1.62039918,
 1.71639094, 1.81999506, 1.93335501, 2.05975875, 2.20481267, 2.37949053,
 2.61309783, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
 0.00000000, 0.00000000, 0.00000000, 0.00000000 ])
varlist['delta']= np.array([0.00946046, 0.01374876, 0.01800757, 0.02223544, 0.026432,   0.03059698,\
 0.03473095, 0.038831,   0.04289777, 0.04693171, 0.05093049, 0.05489528,\
 0.05882575, 0.06272063, 0.06658121, 0.07040614, 0.07419578, 0.07794811,\
 0.08166657, 0.08535033, 0.08900042, 0.09261682, 0.09619758, 0.09974838,\
 0.10326324, 0.10674989, 0.11020544, 0.11363556, 0.11703701, 0.12041313,\
 0.12377785, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,\
 0.00000000, 0.00000000, 0.00000000, 0.00000000])
varlist['R0'] = np.array([10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5,\
 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, \
 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5])
varlist['R1'] = np.array([10.499972,10.499913,10.499826,10.49971 ,10.499565,10.499391,10.49919 ,\
 10.498961,10.498703,10.498419,10.498107,10.497768,10.497403,10.497012,\
 10.496596,10.496154,10.495688,10.495198,10.494685,10.494148,10.493589,\
 10.493009,10.492407,10.491785,10.491142,10.490481,10.489801,10.489103,\
 10.488388,10.487656,10.486909,10.486146,10.485369,10.484579,10.483775,\
 10.482959,10.482131,10.481292,10.480442,10.479583])
varlist['R100'] = np.array([10.499684,10.498987,10.498003,10.496821,10.49557 ,10.494417,10.49353 ,\
 10.493086,10.49324 ,10.494061,10.495527,10.497414,10.499206,10.499996,\
 10.498346,10.492196,10.478827,10.455062,10.41772 ,10.364422,10.294694,\
 10.21081 ,10.118173,10.024561, 9.9389  , 9.870215, 9.828077, 9.82713 ,\
  9.902778,10.161439,10.620902, 9.51475 ,10.393219, 9.826615, 9.592683,\
  9.505773, 9.509459, 9.598973, 9.773525,10.012545])
varlist['R500'] = np.array([10.504281,10.504111,10.501835,10.496888,10.488587,10.476658,10.462013,\
 10.447162,10.436594,10.435952,10.449685,10.475808,10.498375,10.48649 ,\
 10.417483,10.315985,10.249916,10.288231,10.455444,10.583435,10.368276,\
 10.086498,10.103007,10.451723,10.226203, 9.870922,10.350237,10.034389,\
  9.831472, 9.803151,10.467593, 9.478629, 9.697921,10.407131, 9.636414,\
  9.927026, 9.672292,10.277331, 9.50395 ,10.368171])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    print(np.array2string(cal,precision=6, separator=','))
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    perct = 100*sum(abs(act-cal)/div)
    print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




