#!/usr/bin/env python3
import sys, ossys.path.insert(0, '../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.libstell import read_vmec

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)
failtol = 1.0
filename='wout_DSHAPE.nc'
data=read_vmec(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)
else:
    print('EXTENSION: '+filename)
print('==== Scalars ====')
varlist={}
varlist['aspect']=2.9277749505154276
varlist['b0']=0.20584681951892525
varlist['wp']=0.001808732590778445
varlist['betatot']=0.029160921095561083
varlist['volume']=99.45700630158791
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

#print(data['xn'].tolist())

varlist['iotaf']=np.array([[1.0], [0.9947244094488188], [0.9894488188976378], [0.9841732283464566], [0.9788976377952756], [0.9736220472440944], [0.9683464566929134], [0.9630708661417322], [0.9577952755905512], [0.95251968503937], [0.947244094488189], [0.9419685039370078], [0.9366929133858268], [0.9314173228346456], [0.9261417322834646], [0.9208661417322834], [0.9155905511811024], [0.9103149606299212], [0.9050393700787402], [0.899763779527559], [0.894488188976378], [0.8892125984251968], [0.8839370078740157], [0.8786614173228346], [0.8733858267716534], [0.8681102362204725], [0.8628346456692912], [0.8575590551181103], [0.852283464566929], [0.8470078740157481], [0.8417322834645669], [0.8364566929133859], [0.8311811023622047], [0.8259055118110237], [0.8206299212598425], [0.8153543307086615], [0.8100787401574803], [0.8048031496062993], [0.7995275590551181], [0.794251968503937], [0.7889763779527559], [0.7837007874015749], [0.7784251968503937], [0.7731496062992126], [0.7678740157480315], [0.7625984251968504], [0.7573228346456693], [0.7520472440944881], [0.7467716535433071], [0.741496062992126], [0.7362204724409448], [0.7309448818897637], [0.7256692913385827], [0.7203937007874015], [0.7151181102362205], [0.7098425196850393], [0.7045669291338583], [0.6992913385826771], [0.694015748031496], [0.6887401574803149], [0.6834645669291339], [0.6781889763779527], [0.6729133858267716], [0.6676377952755905], [0.6623622047244094], [0.6570866141732283], [0.6518110236220472], [0.6465354330708661], [0.641259842519685], [0.6359842519685039], [0.6307086614173227], [0.6254330708661417], [0.6201574803149605], [0.6148818897637796], [0.6096062992125983], [0.6043307086614174], [0.5990551181102363], [0.593779527559055], [0.588503937007874], [0.583228346456693], [0.5779527559055118], [0.5726771653543308], [0.5674015748031496], [0.5621259842519686], [0.5568503937007874], [0.5515748031496063], [0.5462992125984252], [0.5410236220472441], [0.535748031496063], [0.5304724409448819], [0.5251968503937008], [0.5199212598425197], [0.5146456692913386], [0.5093700787401574], [0.5040944881889764], [0.49881889763779524], [0.49354330708661415], [0.48826771653543305], [0.48299212598425195], [0.4777165354330708], [0.47244094488188965], [0.4671653543307086], [0.46188976377952756], [0.4566141732283464], [0.4513385826771653], [0.44606299212598427], [0.4407874015748031], [0.435511811023622], [0.4302362204724409], [0.42496062992125977], [0.41968503937007867], [0.4144094488188976], [0.4091338582677165], [0.4038582677165354], [0.3985826771653543], [0.3933070866141732], [0.3880314960629921], [0.382755905511811], [0.37748031496062984], [0.37220472440944874], [0.3669291338582677], [0.36165354330708654], [0.35637795275590545], [0.3511023622047244], [0.34582677165354325], [0.3405511811023621], [0.335275590551181], [0.3299999999999999]])
varlist['presf']=np.array([[1599.9255998511992], [1574.9271498542998], [1550.0279000558005], [1525.3270506541014], [1500.8246016492035], [1476.5205530411063], [1452.4149048298098], [1428.5076570153144], [1404.7988095976193], [1381.2883625767251], [1357.976315952632], [1334.8626697253396], [1311.947423894848], [1289.230578461157], [1266.712133424267], [1244.3920887841778], [1222.270444540889], [1200.3472006944016], [1178.6223572447147], [1157.0959141918286], [1135.767871535743], [1114.6382292764586], [1093.706987413975], [1072.974145948292], [1052.4397048794099], [1032.1036642073286], [1011.966023932048], [992.0267840535682], [972.2859445718892], [952.7435054870112], [933.3994667989338], [914.2538285076571], [895.3065906131812], [876.5577531155062], [858.0073160146322], [839.6552793105587], [821.5016430032861], [803.5464070928143], [785.7895715791431], [768.2311364622728], [750.8711017422035], [733.709467418935], [716.746233492467], [699.9813999628], [683.4149668299337], [667.0469340938682], [650.8773017546036], [634.9060698121397], [619.1332382664766], [603.5588071176144], [588.1827763655527], [573.0051460102919], [558.025916051832], [543.245086490173], [528.6626573253147], [514.2786285572571], [500.0930001860005], [486.10577221154455], [472.3169446338894], [458.72651745303494], [445.3344906689814], [432.1408642817287], [419.14563829127667], [406.3488126976255], [393.7503875007751], [381.35036270072555], [369.14873829747677], [357.1455142910287], [345.3406906813814], [333.7342674685348], [322.32624465248915], [311.11662223324436], [300.1054002108003], [289.2925785851572], [278.67815735631467], [268.262136524273], [258.04451608903213], [248.02529605059206], [238.20447640895287], [228.5820571641144], [219.15803831607667], [209.93241986483977], [200.9052018104037], [192.0763841527684], [183.44596689193386], [175.01395002790017], [166.78033356066726], [158.74511749023512], [150.9083018166037], [143.26988653977307], [135.82987165974333], [128.58825717651453], [121.54504309008632], [114.70022940045882], [108.05381610763219], [101.60580321160637], [95.3561907123814], [89.30497860995717], [83.45216690433378], [77.79775559551112], [72.34174468348922], [67.08413416826816], [62.024924049847925], [57.16411432822852], [52.501705003409874], [48.03769607539197], [43.772087544174994], [39.70487940975875], [35.83607167214326], [32.165664331328614], [28.69365738731471], [25.42005084010173], [22.34484468968949], [19.468038936077917], [16.789633579267264], [14.309628619257353], [12.028024056048194], [9.944819889639957], [8.060016120032378], [6.373612747225633], [4.885609771219723], [3.5960071920144725], [2.5048050096101453], [1.6120032240066529], [0.9176018352038184], [0.42160084320164165], [0.12400024800038878], [-0.0744001488002688]])
varlist['jcurv']=np.array([[-70186.81446507735], [-69697.37255821451], [-69207.93065135166], [-68717.13777851399], [-68225.05941168209], [-67731.72226784001], [-67236.55259976214], [-66739.57821125173], [-66240.9012387713], [-65740.64385938004], [-65238.888329542875], [-64735.66484651119], [-64230.97906265901], [-63724.82225711096], [-63217.17675620758], [-62708.023575785744], [-62197.34800586784], [-61685.137760818936], [-61171.38135437527], [-60656.06711422136], [-60139.185263166284], [-59620.72710603025], [-59100.68545203823], [-58579.052670570876], [-58055.82187728137], [-57530.98582499379], [-57004.53883492611], [-56476.47531731125], [-55946.79133793153], [-55415.48267485862], [-54882.546645303264], [-54347.98019289876], [-53811.781936498555], [-53273.95029212484], [-52734.4854811876], [-52193.3874906001], [-51650.6582506762], [-51106.29943505572], [-50560.31465801258], [-50012.707394440826], [-49463.48307183392], [-48912.64715083216], [-48360.20705750745], [-47806.1703140187], [-47250.5465515709], [-46693.345586920695], [-46134.57931240486], [-45574.26007983419], [-45012.40225141392], [-44449.02084986638], [-43884.13292252268], [-43317.75646867], [-42749.91141582445], [-42180.618977584076], [-41609.90222798641], [-41037.78587241848], [-40464.296439675636], [-39889.462381213416], [-39313.31394353283], [-38735.88378356354], [-38157.20617508884], [-37577.318237185435], [-36996.25876309898], [-36414.069733842334], [-35830.794774538146], [-35246.481192422594], [-34661.1780031648], [-34074.93825495005], [-33487.81681522285], [-32899.87297149485], [-32311.168198712472], [-31721.768442928973], [-31131.742291808034], [-30541.163026048347], [-29950.107391634225], [-29358.656161328385], [-28766.894867074967], [-28174.912874048598], [-27582.80571641727], [-26990.671294504788], [-26398.615994441498], [-25806.747666175033], [-25215.18503808686], [-24624.045400646006], [-24033.460807475738], [-23443.559794801055], [-22854.48980292641], [-22266.389708016668], [-21679.423636515723], [-21093.742726704942], [-20509.530109307678], [-19926.947737701284], [-19346.199201658776], [-18767.459245727532], [-18190.95530088192], [-17616.873921822713], [-17045.467298231906], [-16476.933366153637], [-15911.5536592725], [-15349.538325994872], [-14791.197880078542], [-14236.752005129485], [-13686.547795439825], [-13140.816269691337], [-12599.93520155246], [-12064.140421513768], [-11533.857654924957], [-11009.333203188093], [-10491.014330064996], [-9979.149099508228], [-9474.24880400056], [-8976.5631781677], [-8486.608156877644], [-8004.650888113412], [-7531.283402178399], [-7066.733953233381], [-6611.544144316689], [-6166.019755690035], [-5730.897168334461], [-5306.153416968977], [-4892.564226395024], [-4489.809999665997], [-4101.315847822883], [-3721.7794298034705], [-3358.917161372792], [-3005.989286285825], [-2665.594791527809], [-2325.2002967697927]])
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    perct = 100*sum(abs(act-cal)/act)
    print('  '+temp+': '+str(cal[0])+'   '+str(act[0])+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('==== Arrays =====')
varlist={}
varlist['rmnc']=np.array([0.3881599981989999, 0.543839182013267,0.7531583987627359])
varlist['zmns']=np.array([0.498612043832715,0.7068949306570961,1.0066169239309073])
varlist['lmns']=np.array([-0.24066933660557893,-0.34209941038013136,-0.48880011719128896])
varlist['bmnc']=np.array([-0.02166836797715306,-0.030645003303238893,-0.04350649740090247])
for temp in varlist:
    mn  = 1; k   = 16
    act = varlist[temp][0]
    cal = data[temp][k,mn]
    s   = str(k)
    m   = str(int(data['xm'][mn]))
    n   = str(int(data['xn'][mn]))
    perct = 100*abs((act-cal)/act)
    print('  '+temp+'['+s+','+m+','+n+']: '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
    mn  = 1; k   = 32
    act = varlist[temp][1]
    cal = data[temp][k,mn]
    s   = str(k)
    m   = str(int(data['xm'][mn]))
    n   = str(int(data['xn'][mn]))
    perct = 100*abs((act-cal)/act)
    print('  '+temp+'['+s+','+m+','+n+']: '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
    mn  = 1; k   = 64
    act = varlist[temp][2]
    cal = data[temp][k,mn]
    s   = str(k)
    m   = str(int(data['xm'][mn]))
    n   = str(int(data['xn'][mn]))
    perct = 100*abs((act-cal)/act)
    print('  '+temp+'['+s+','+m+','+n+']: '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)




