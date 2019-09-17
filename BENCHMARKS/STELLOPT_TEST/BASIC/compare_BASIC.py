#!/usr/bin/env python3
import sys, ossys.path.insert(0, '../../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.stellopt import read_stellopt

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

failtol = 1.0
filename='stellopt.BASIC'
data=read_stellopt(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)
else:
    print('EXTENSION: '+filename)
print('==== Scalars ====')
varlist={}
varlist['PHIEDGE_equil']=0.514386
varlist['CURTOR_equil']=-1.742957864085E+005
varlist['RBTOR_equil']=2.312873858830
varlist['R0_equil']=1.577564814144
varlist['VOLUME_equil']=2.9787172145375
varlist['BETA_equil']=4.247495078898E-002
varlist['WP_equil']=1.925493826145E+005
varlist['ASPECT_equil']=4.365254725961
lfail = 0;
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    perct = 100*(abs(act-cal)/act)
    print('  '+temp+': '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('==== Vectors ====')
varlist={}

varlist['NELINE_equil']=np.array([3.234863326891e+19, 4.207356830195e+19, 4.942220121813e+19, 4.971516997227e+19, 4.951985746951e+19, 4.989420643313e+19, 4.960123767899e+19, 4.989420643313e+19, 4.950358142762e+19, 4.969889393037e+19, 4.979655018175e+19, 4.989420643313e+19, 4.989420643313e+19, 4.979655018175e+19, 4.961751372089e+19, 4.932454496676e+19, 4.969889393037e+19, 4.991048247503e+19, 4.9129232464e+19, 4.991048247502e+19, 4.950358142762e+19, 4.969889393037e+19, 4.940592517624e+19, 4.951985746951e+19, 4.991048247502e+19, 4.922688871538e+19, 4.979655018175e+19, 0.0, 0.0, 0.0])
varlist['FARADAY_equil']=np.array([7.139215399803e+18, 9.052056823689e+18, 1.054155604138e+19, 1.056794570172e+19, 1.0456362497e+19, 1.017734323199e+19, 1.009362064842e+19, 9.824986289389e+18, 9.526748229683e+18, 9.378857918157e+18, 9.086778408452e+18, 8.866539545873e+18, 8.605215214709e+18, 8.242644295106e+18, 8.026901997195e+18, 7.694140815456e+18, 7.343042310993e+18, 7.009000035417e+18, 6.633350729282e+18, 6.264483556361e+18, 5.845752589955e+18, 5.46653316302e+18, 4.899797567804e+18, 4.383240651682e+18, 3.680359872178e+18, 2.794000986661e+18, 1.384466723698e+18, 0.0, 0.0, 0.0])
varlist['PRESS_equil']=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.003560065575185, 0.1031741697915, 0.2733332450359, 0.4405628979824, 0.5771797184444, 0.6823551168574, 0.7624073229203, 0.8236392002169, 0.8707098074245, 0.9067670716376, 0.9339953408312, 0.9540832675538, 0.9684800034966, 0.9784838824074, 0.9852344842785, 0.9896804531468, 0.9925510894581, 0.9943719214941, 0.9954841942988, 0.9961029217696, 0.9963654916206, 0.9963591833772, 0.9961593569576, 0.9958339467206, 0.9954438835409, 0.9950402655701, 0.9946756485326, 0.9943997563471, 0.9942512135293, 0.9942512140801, 0.9943997579482, 0.9946756510376, 0.9950402687614, 0.9954438871657, 0.9958339505838, 0.9961593608753, 0.9963591871799, 0.9963654951787, 0.9961029249774, 0.9954841970368, 0.9943719236116, 0.9925510905556, 0.9896804525408, 0.9852344808764, 0.9784838746256, 0.9684799892414, 0.9540832442439, 0.9339953053966, 0.9067670204766, 0.8707097360734, 0.8236391026139, 0.7624071901648, 0.6823549358362, 0.5771794720043, 0.4405625711729, 0.2733328516103, 0.1031738183249, 0.003559983065329, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
varlist['NE_equil']=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
varlist['TE_equil']=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.72535141048, 471.9346609297, 1251.161798764, 2016.861833142, 2642.328553287, 3123.801034882, 3490.246584935, 3770.537216928, 3986.010840727, 4151.075868028, 4275.729194636, 4367.694630125, 4433.605235154, 4479.400532239, 4510.29898452, 4530.640380848, 4543.776921854, 4552.099612378, 4557.190921618, 4560.035324749, 4561.236074888, 4561.207243668, 4560.29358245, 4558.80212399, 4557.00527658, 4555.153708361, 4553.487281259, 4552.226843861, 4551.547617573, 4551.547620092, 4552.226851179, 4553.487292701, 4555.153722962, 4557.005293266, 4558.802141736, 4560.293600374, 4561.207261049, 4561.23609115, 4560.035339433, 4557.190934234, 4552.099622058, 4543.776926867, 4530.640378074, 4510.298968946, 4479.400496614, 4433.605169893, 4367.694523391, 4275.729032416, 4151.075633815, 3986.010514103, 3770.536770139, 3490.245977189, 3123.800206228, 2642.327425121, 2016.860336457, 1251.15999803, 471.9330526949, 15.72497414425, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
varlist['TI_equil']=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 419.3295522232, 1399.265659004, 2318.863506893, 3010.402649893, 3500.700945034, 3850.597354456, 4103.298319454, 4283.677295642, 4407.319425897, 4486.677753339, 4533.381483685, 4558.298275196, 4570.466352163, 4575.803433794, 4576.953953376, 4574.166220348, 4565.04416692, 4541.399347934, 4484.198416614, 4363.27065221, 4146.049025335, 3803.274721578, 3297.424723466, 2570.704607261, 1583.83503038, 480.4338273235, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
varlist['BALLOON_grate']=np.array([-0.0399052236745, -0.04518164162041, -0.04947750299293, -0.04971083018054, -0.05074735892034, -0.05734150402709, -0.06015735850861, -0.06090566811739, -0.06153402446057, -0.06276648035781, -0.06350834493599, -0.06358017707733, -0.06444598930457, -0.06866249806963, -0.07489275402092, -0.08117312954453, -0.0861868793103, -0.09018186356229, -0.09306345912829, -0.09518753535743, -0.09651640408053, -0.09723828136738, -0.09734683810026, -0.09696629216251, -0.09622157938765, -0.09536946729703, -0.09483969201022, -0.09523124002104, -0.09709229043953, -0.1005133158428, -0.1048982645875, -0.1095641363109, -0.113872228696, -0.1176107275844, -0.1206010750882, -0.1229150547474, -0.124484269068, -0.1254064183285, -0.1256202823846, -0.1252147741061, -0.1241314706584, -0.1224541243388, -0.1201472542902, -0.1173096555908, -0.1139874689766, -0.1104500109516, -0.1071520857077, -0.1049573902257, -0.1048809461717, -0.107615407385, -0.1130019491296, -0.1200346076679, -0.1274577349067, -0.1343054407947, -0.1401069712201, -0.1446471288074, -0.148032228357, -0.1502651638187, -0.1516290629516, -0.1521303609038, -0.152112512322, -0.1514716723326, -0.1505717755191, -0.1492286912724, -0.1478275417753, -0.1459996078021, -0.1441455282454, -0.1420397585565, -0.1397478834653, -0.1372387769706, -0.1348156464585, -0.1319618852003, -0.1297473360322, -0.1274352361651, -0.1268689225235, -0.1253819594997, -0.1247349886582, -0.1176968585891, -0.1106046362363, -0.08907185574398, -0.07101455369933, 0.026928735503, 0.06124562207265, 0.09033338664691, 0.08985874270061, 0.08501649973644, 0.04413508386684, -0.05680148431873, -0.07107962259673, 0.04761966498712, 0.07757838930801, -0.08986066241933, -0.1363185063141, -0.1480696220691, -0.1508076730747, -0.1489653063131, -0.1443836841124, 0.0])
varlist['BOOTSTRAP_equil']=np.array([15167.3233769, 18501.21754991, 24198.02798831, 29256.24018519, 31846.00594443, 20402.6861092, 123054.8940804, 87516.58887283, 90167.27963477, 96666.02737257, 105872.6900355, 114200.4455997, 125481.3309881, 136315.6343484, 145948.802015, 165810.5598879, 172263.0287623, 183355.1654364, 192387.6580427, 216825.2881505, 224590.4328844, 236393.6285925, 248683.9542036, 261566.5202209, 272901.4082687, 284612.2405254, 294218.294021, 300934.2874703, 321790.7915556, 330871.8788079, 338792.508678, 347458.7446072, 357107.5076378, 369484.0086015, 377114.8608438, 376393.46389, 397977.4176806, 400242.7577104, 442188.236557, 426086.7568332, 452142.2674224, 453617.4971625, 454063.14352, 478105.279068, 430076.0954915, 401979.2655309, 371379.297407, 363535.3432868, 304100.6980056, 405152.059098, 308347.9557527, 241084.0452573, 261893.4163625, 160376.7144049, 118800.0975831, 31296.01670631, 94793.84970213, 304688.0736922, 703386.9882114, 1707738.470246, 6952039.497356, 3932083.640475, 2053270.967182, 1471520.281892, 1196202.361154, 1024141.992902, 917540.7037212, 846111.8893132, 714997.4895081, 704187.1188931, 591587.3007793, 673049.7947525, 574466.7740582, 522471.4891911, 392249.8106689, 723543.3333808, 546439.4031751, 470307.9802537, 421811.0307938, 408608.3234336, 382927.6361731, 339154.1847212, 324220.7765603, 288983.2662726, 261852.7384881, 233963.5427849, 206567.245802, 181422.8249552, 158247.1496434, 137226.6951551, 118073.624366, 100134.877888, 82698.17961471, 64874.36949482, 46058.69134777, 1.81904195909e+32, 4.798274102878e+35, 0.0002985394726398])
varlist['NEO_equil']=np.array([1.634137105073e-07, 9.809754629811e-08, 1.400181030631e-07, 2.232633990728e-07, 3.348331508397e-07, 4.771001541935e-07, 6.544552711931e-07, 8.463282409109e-07, 1.106276538918e-06, 1.401126741018e-06, 1.739873387997e-06, 2.070469999115e-06, 2.499767064814e-06, 2.958590913581e-06, 3.524093624128e-06, 4.115223860603e-06, 4.789509674426e-06, 5.553866057218e-06, 6.171442557097e-06, 7.251638124726e-06, 8.165732167825e-06, 9.237593400606e-06, 1.042747606882e-05, 1.171033248935e-05, 1.300905871199e-05, 1.422415877213e-05, 1.615751531158e-05, 1.787123982441e-05, 1.97005177205e-05, 2.186407889346e-05, 2.461144036412e-05, 2.652980343614e-05, 2.954979814647e-05, 3.248697972839e-05, 3.548133352086e-05, 3.907233210005e-05, 4.293630150512e-05, 4.675452792637e-05, 5.076569322065e-05, 5.551921841857e-05, 6.012558232933e-05, 6.571286283634e-05, 7.125794800456e-05, 7.705449733924e-05, 8.390622745592e-05, 9.026687647197e-05, 9.740492986519e-05, 0.0001050359961178, 0.0001134283127953, 0.000122632880009, 0.0001312556856943, 0.0001413763963709, 0.0001514085372828, 0.0001638336544956, 0.0001756014619506, 0.0001873338916894, 0.0001998935172734, 0.0002131031268064, 0.0002268680160154, 0.0002423056672998, 0.0002580251395367, 0.0002753175413953, 0.0002917057159331, 0.0003113101845538, 0.0003287394354777, 0.0003467046644809, 0.0003684276661603, 0.0003898284807443, 0.0004126209410499, 0.0004361993263121, 0.0004605802501365, 0.0004871364700153, 0.0005110434312719, 0.0005435483345191, 0.0005736738665381, 0.0006069066058251, 0.0006442646223093, 0.0006823220079323, 0.0007211824477164, 0.0007622353812488, 0.0008113963238251, 0.0008570113619511, 0.0009085047033222, 0.0009695787350655, 0.001025351836078, 0.001092259760887, 0.001163982649923, 0.00124113111335, 0.001309488486633, 0.001393117461776, 0.001494539150766, 0.001587732897119, 0.00167897461567, 0.001800374052995, 0.001892534830587, 0.002010384130923, 0.002118903614155, 0.002263698301993])
varlist['TXPORT_equil']=np.array([0.09589212091657, 0.09667090971904, 0.09587599783759, 0.09536676885799, 0.0962629335895, 0.09797820318666, 0.09970560774425, 0.1007944843387, 0.09998994522488, 0.09477010520523, 0.09064250460416, 0.08752356857393, 0.0849144207205, 0.08499795558171, 0.09148810219251])

#print(data['HELICITY_FULL_equil'][25].tolist())
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    perct = 100*sum(abs(act-cal)/max(abs(act)))
    print('  '+temp+': '+str(max(cal))+'   '+str(max(act))+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




