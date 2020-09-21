Tutorial: VMEC Fixed Boundary Run
=================================

This tutorial will walk the user through running VMEC with a fixed
boundary condition. For this example the National Compact Stellarator
Experiment (NCSX) configuration will be used. This machine is
stellarator symmetric with a periodicity of three.

------------------------------------------------------------------------

1.  \> \_\_**Edit the input namelist text file.**\_\_ \> ﻿The input
    namelist (<file:input.ncsx_c09r00_fixed>) controls the execution of
    the VMEC code. The suffix of the input file will be appended to each
    of the output files as we will see after execution. The Fourier
    coefficient in this file have been generated through an optimization
    routine. In general, more simple initial conditions will suffice for
    the axis position and outer most flux surface.

<!-- -->

    #!fortran
    &INDATA
    !----- Runtime Parameters -----
      DELT =   9.00E-01
      NITER = 2500
      NSTEP =  200
      TCON0 =   2.00E+00
      NS_ARRAY =  9   28   49
      FTOL_ARRAY =  1.000000E-06  1.000000E-08  1.000000E-11
    !----- Grid Parameters -----
      LASYM = F
      NFP =    3
      MPOL =   11
      NTOR =    6
      NZETA  =   32
      PHIEDGE =   4.97070205837336E-01
    !----- Free Boundary Parameters -----
      LFREEB = F
      MGRID_FILE = 'none'
      EXTCUR = 0
      NVACSKIP =    0
    !----- Pressure Parameters -----
      GAMMA =   0.000000E+00
      BLOAT =   1.000000E+00
      SPRES_PED =   1.00000000000000E+00
      AM =   6.85517649352426E+04 -5.12027745123057E+03 -3.61510451745464E+04 -4.74263014113066E+05
      1.78878195473870E+06 -3.21513828868170E+06  2.69041023837233E+06 -8.17049854168367E+05
      0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
    !----- Current/Iota Parameters -----
      CURTOR =  -1.78606250000000E+05
      NCURR =    1
      AI =   0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
      0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
      0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
      AC =   8.18395699999999E+03  1.43603560000000E+06 -1.07407140000000E+07  7.44389200000000E+07
     -3.22215650000000E+08  8.81050800000000E+08 -1.49389660000000E+09  1.52746800000000E+09
     -8.67901590000000E+08  2.10351200000000E+08  0.00000000000000E+00
    !----- Axis Parameters -----
      RAXIS =   1.49569454253276E+00  1.05806400912320E-01  7.21255454715878E-03 -3.87402652289249E-04
     -2.02425864534069E-04 -1.62602353744308E-04 -8.89569831063077E-06
      ZAXIS =   0.00000000000000E+00 -5.19492027001782E-02 -3.18814224021375E-03  2.26199929262002E-04
      1.28803681387330E-04  1.11266150452637E-06  1.13732703961869E-05
    !----- Boundary Parameters -----
      RBC(0,0) =   1.40941668895656E+00     ZBS(0,0) =   0.00000000000000E+00
      RBC(1,0) =   2.79226697269945E-02     ZBS(1,0) =  -1.92433268059631E-02
      RBC(2,0) =  -1.54739398509667E-03     ZBS(2,0) =   1.11459511078088E-02
      RBC(3,0) =   2.90733840104882E-03     ZBS(3,0) =  -3.97869471888770E-03
      RBC(4,0) =  -8.91322016448873E-04     ZBS(4,0) =   1.34394804673514E-03
      RBC(5,0) =  -7.81997770407839E-05     ZBS(5,0) =  -1.57143910159387E-04
      RBC(6,0) =   1.06129711928351E-04     ZBS(6,0) =   9.58024291307491E-05
      RBC(-6,1) =  2.48228899767757E-05     ZBS(-6,1) = -2.28386224209054E-05
      RBC(-5,1) =  8.23567667077671E-05     ZBS(-5,1) =  3.30176003890210E-04
      RBC(-4,1) = -7.20220898033597E-04     ZBS(-4,1) =  1.28038328362904E-04
      RBC(-3,1) =  2.76250777733235E-03     ZBS(-3,1) =  3.43199911886317E-04
      RBC(-2,1) = -1.24883373588382E-02     ZBS(-2,1) =  6.12174680232785E-04
      RBC(-1,1) =  1.52272804511910E-02     ZBS(-1,1) = -2.70066914159594E-02
      RBC(0,1) =   2.89195233044040E-01     ZBS(0,1) =   4.50462554508443E-01
      RBC(1,1) =  -1.17988850341728E-01     ZBS(1,1) =   1.93490230971634E-01
      RBC(2,1) =  -3.84923299492945E-03     ZBS(2,1) =   5.72865331625290E-03
      RBC(3,1) =  -1.44452305429529E-03     ZBS(3,1) =   2.19788951889214E-03
      RBC(4,1) =  -2.11622985211109E-04     ZBS(4,1) =   1.31883972780290E-03
      RBC(5,1) =   1.79091719677667E-04     ZBS(5,1) =  -5.63363462408534E-04
      RBC(6,1) =   1.31982402741742E-04     ZBS(6,1) =  -9.31801467009349E-05
      RBC(-6,2) = -2.40882614870476E-05     ZBS(-6,2) = -3.95416405717970E-05
      RBC(-5,2) = -4.92449386382591E-05     ZBS(-5,2) = -3.25048356502217E-06
      RBC(-4,2) =  1.50530476034115E-04     ZBS(-4,2) =  4.61522421935086E-05
      RBC(-3,2) = -1.23084235126550E-03     ZBS(-3,2) = -3.40868203306282E-04
      RBC(-2,2) =  2.01350576071929E-04     ZBS(-2,2) = -4.19781517712033E-03
      RBC(-1,2) =  2.36777003797179E-03     ZBS(-1,2) =  1.98753868216412E-02
      RBC(0,2) =   5.73443941583452E-02     ZBS(0,2) =   4.81527027892127E-03
      RBC(1,2) =   6.89385874058265E-02     ZBS(1,2) =  -9.28353553039424E-03
      RBC(2,2) =   4.71996849673782E-02     ZBS(2,2) =  -2.04292782322197E-02
      RBC(3,2) =  -5.50889052720066E-04     ZBS(3,2) =   8.81593501270446E-04
      RBC(4,2) =   4.24491391207156E-04     ZBS(4,2) =  -6.08871281835245E-04
      RBC(5,2) =  -2.07538883155595E-04     ZBS(5,2) =  -3.88708113241096E-04
      RBC(6,2) =  -1.62304038006678E-04     ZBS(6,2) =   1.72340342752605E-04
      RBC(-6,3) = -1.01105699684233E-04     ZBS(-6,3) = -6.16215454248342E-05
      RBC(-5,3) =  5.15925605980565E-05     ZBS(-5,3) =  1.23419431936950E-04
      RBC(-4,3) = -3.79290487874111E-05     ZBS(-4,3) =  3.98637008165582E-06
      RBC(-3,3) = -2.96154201246223E-04     ZBS(-3,3) = -7.01248486620889E-04
      RBC(-2,3) =  1.27628943631957E-03     ZBS(-2,3) =  3.19332333533202E-03
      RBC(-1,3) =  3.12803506573940E-03     ZBS(-1,3) = -8.24657727838880E-03
      RBC(0,3) =  -1.34574092972690E-02     ZBS(0,3) =   5.05936199755365E-03
      RBC(1,3) =  -8.02339287294677E-03     ZBS(1,3) =  -3.90421394288867E-03
      RBC(2,3) =  -1.68510947837154E-02     ZBS(2,3) =   3.75441853342170E-03
      RBC(3,3) =  -8.00581733372124E-03     ZBS(3,3) =   6.00542774606014E-03
      RBC(4,3) =   1.80667899211621E-03     ZBS(4,3) =  -4.16787432635077E-04
      RBC(5,3) =   3.10773970094350E-05     ZBS(5,3) =   5.44335921432213E-05
      RBC(6,3) =   8.32496816115997E-05     ZBS(6,3) =  -4.15830451164888E-05
      RBC(-6,4) = -1.19874891436340E-05     ZBS(-6,4) =  1.56845408711308E-05
      RBC(-5,4) =  1.22793444338155E-04     ZBS(-5,4) = -3.97576733690054E-05
      RBC(-4,4) = -1.30945484439682E-04     ZBS(-4,4) = -7.22429623460448E-05
      RBC(-3,4) = -1.21368603604647E-04     ZBS(-3,4) =  3.52928331257216E-04
      RBC(-2,4) =  1.00352526472782E-03     ZBS(-2,4) = -1.23710282249961E-04
      RBC(-1,4) = -1.73680844498789E-03     ZBS(-1,4) = -1.50689928334813E-03
      RBC(0,4) =   1.80149787198970E-03     ZBS(0,4) =   1.56109492686192E-03
      RBC(1,4) =   3.82771889154294E-03     ZBS(1,4) =   3.80910842862487E-03
      RBC(2,4) =   5.43835034437129E-03     ZBS(2,4) =   2.06275075117804E-03
      RBC(3,4) =   8.39729828422411E-04     ZBS(3,4) =  -1.54779126563731E-03
      RBC(4,4) =   6.74263596810560E-04     ZBS(4,4) =  -1.33149943553452E-03
      RBC(5,4) =  -6.98647584180715E-04     ZBS(5,4) =   3.81307095116973E-04
      RBC(6,4) =   8.77670652920776E-05     ZBS(6,4) =  -1.40433963574141E-05
      RBC(-6,5) =  6.78635213884316E-06     ZBS(-6,5) = -1.22283666932084E-05
      RBC(-5,5) =  3.87846546342867E-05     ZBS(-5,5) =  4.64829761643373E-05
      RBC(-4,5) = -3.78300368387435E-05     ZBS(-4,5) = -7.03801581329045E-05
      RBC(-3,5) = -1.21743926248229E-05     ZBS(-3,5) =  1.85735151533626E-04
      RBC(-2,5) = -2.68229697014545E-04     ZBS(-2,5) = -9.33216243296025E-04
      RBC(-1,5) =  1.19567316567517E-03     ZBS(-1,5) =  2.12648562837673E-03
      RBC(0,5) =  -7.12579133390599E-04     ZBS(0,5) =  -1.97890515574565E-03
      RBC(1,5) =   8.81127157923892E-04     ZBS(1,5) =   2.71321673191593E-03
      RBC(2,5) =   9.67210453659238E-04     ZBS(2,5) =   8.74618447862515E-04
      RBC(3,5) =   2.11794179698155E-04     ZBS(3,5) =   8.43817701627930E-04
      RBC(4,5) =   1.29403911922840E-03     ZBS(4,5) =   6.51808476607835E-04
      RBC(5,5) =  -1.30477683585083E-04     ZBS(5,5) =   1.01349326961770E-04
      RBC(6,5) =   1.86680624010370E-04     ZBS(6,5) =  -2.13838628730300E-04
      RBC(-6,6) = -4.08213549686361E-05     ZBS(-6,6) = -7.53394429655583E-06
      RBC(-5,6) =  7.11305157811999E-05     ZBS(-5,6) =  2.54876062250879E-05
      RBC(-4,6) = -1.33727065581923E-04     ZBS(-4,6) = -1.70180862196520E-05
      RBC(-3,6) =  1.65191943182183E-06     ZBS(-3,6) = -1.31350577800873E-04
      RBC(-2,6) =  2.19460449719541E-04     ZBS(-2,6) =  4.38914760402648E-04
      RBC(-1,6) =  4.68618562605432E-04     ZBS(-1,6) = -4.44537659614533E-04
      RBC(0,6) =  -8.51896573200937E-04     ZBS(0,6) =   7.36122964253313E-04
      RBC(1,6) =  -5.26623264534578E-05     ZBS(1,6) =  -1.12352425125337E-03
      RBC(2,6) =  -1.31954654361710E-04     ZBS(2,6) =  -2.22905186553194E-03
      RBC(3,6) =  -8.91482312658694E-04     ZBS(3,6) =  -2.11193996461398E-03
      RBC(4,6) =  -3.89733094884781E-04     ZBS(4,6) =  -3.44184359663702E-04
      RBC(5,6) =  -2.74329775462215E-04     ZBS(5,6) =  -5.06914660659672E-05
      RBC(6,6) =   2.47385092660320E-04     ZBS(6,6) =   3.74971583066409E-05
      RBC(-6,7) =  9.61516193308531E-06     ZBS(-6,7) = -3.66121037894761E-06
      RBC(-5,7) = -2.51122684780459E-05     ZBS(-5,7) =  3.72828134065079E-05
      RBC(-4,7) =  4.44568599556351E-05     ZBS(-4,7) = -8.74488353626824E-05
      RBC(-3,7) = -1.42433799354752E-04     ZBS(-3,7) =  1.48694485468843E-04
      RBC(-2,7) =  4.85802385952487E-04     ZBS(-2,7) = -2.27519962800893E-04
      RBC(-1,7) = -9.00652688032426E-04     ZBS(-1,7) =  4.16601324903870E-04
      RBC(0,7) =   9.59457670863182E-04     ZBS(0,7) =  -3.25818663499641E-04
      RBC(1,7) =  -3.37159659594826E-04     ZBS(1,7) =  -2.34240245561361E-04
      RBC(2,7) =  -4.64969900861713E-04     ZBS(2,7) =   4.87821281121050E-04
      RBC(3,7) =  -4.09185322970312E-04     ZBS(3,7) =   8.50140634573578E-04
      RBC(4,7) =   5.32088748759921E-05     ZBS(4,7) =   5.93528572346752E-04
      RBC(5,7) =  -3.21692982976907E-04     ZBS(5,7) =  -2.54775193277671E-04
      RBC(6,7) =  -4.82403633897412E-05     ZBS(6,7) =   1.41947169759239E-05
      RBC(-6,8) = -2.23522770283961E-05     ZBS(-6,8) = -4.00911971000495E-06
      RBC(-5,8) =  3.95696912099304E-05     ZBS(-5,8) =  1.34684147523625E-05
      RBC(-4,8) = -6.50775924544567E-05     ZBS(-4,8) = -2.94168940555405E-05
      RBC(-3,8) =  1.71610112932980E-04     ZBS(-3,8) =  2.17875987311858E-05
      RBC(-2,8) = -3.45412623614909E-04     ZBS(-2,8) = -7.26482153663716E-05
      RBC(-1,8) =  5.61089095467387E-04     ZBS(-1,8) =  2.51145295676537E-04
      RBC(0,8) =  -5.84359101746051E-04     ZBS(0,8) =  -5.42465826224607E-04
      RBC(1,8) =  -6.16860761080513E-05     ZBS(1,8) =   3.93697603313273E-04
      RBC(2,8) =   5.99275780897287E-04     ZBS(2,8) =   3.30798770955874E-04
      RBC(3,8) =   5.68520162541870E-04     ZBS(3,8) =   5.47788467933391E-04
      RBC(4,8) =   4.47404034542356E-04     ZBS(4,8) =   2.43547539548605E-04
      RBC(5,8) =   2.76704814165950E-04     ZBS(5,8) =   9.15194583315619E-05
      RBC(6,8) =   2.97621090888441E-04     ZBS(6,8) =   1.65605427353701E-04
      RBC(-6,9) = -3.78145897931544E-06     ZBS(-6,9) = -1.85759364750771E-06
      RBC(-5,9) = -1.57985677623482E-05     ZBS(-5,9) =  1.06970045147331E-05
      RBC(-4,9) =  7.91641381274532E-05     ZBS(-4,9) = -1.16252171939772E-05
      RBC(-3,9) = -1.97587428419659E-04     ZBS(-3,9) = -3.08457797690412E-05
      RBC(-2,9) =  3.95855751452672E-04     ZBS(-2,9) =  2.03418980231168E-04
      RBC(-1,9) = -5.41153103221438E-04     ZBS(-1,9) = -3.99552958537408E-04
      RBC(0,9) =   4.98714381092541E-04     ZBS(0,9) =   4.32916759100965E-04
      RBC(1,9) =  -8.06048953531492E-05     ZBS(1,9) =  -1.84722027458208E-04
      RBC(2,9) =  -8.67990109801738E-05     ZBS(2,9) =   2.52568631885491E-04
      RBC(3,9) =   4.35340840113358E-04     ZBS(3,9) =  -3.50159782847442E-04
      RBC(4,9) =   2.33585243788111E-04     ZBS(4,9) =  -7.06133299107118E-04
      RBC(5,9) =  -7.69581174305243E-06     ZBS(5,9) =  -3.79072907220561E-04
      RBC(6,9) =  -2.85256407945938E-05     ZBS(6,9) =  -4.49599333610498E-05
      RBC(-6,10) =  1.20206720198758E-05     ZBS(-6,10) =  4.73580005255806E-06
      RBC(-5,10) =  7.02670536357846E-06     ZBS(-5,10) = -6.99664911015022E-06
      RBC(-4,10) = -2.76926398374910E-05     ZBS(-4,10) = -9.18014408856618E-06
      RBC(-3,10) =  5.20223745639364E-05     ZBS(-3,10) =  6.80574180383753E-05
      RBC(-2,10) = -7.88310431749746E-05     ZBS(-2,10) = -1.06370673487973E-04
      RBC(-1,10) =  1.21755712542490E-05     ZBS(-1,10) =  1.22161894513591E-04
      RBC(0,10) =  3.22193519645521E-05     ZBS(0,10) = -6.04052049877600E-05
      RBC(1,10) = -1.08453911913102E-06     ZBS(1,10) =  8.60890353665843E-05
      RBC(2,10) =  1.04051545504927E-04     ZBS(2,10) = -2.17661420286656E-04
      RBC(3,10) = -5.21965328013036E-04     ZBS(3,10) = -2.67111216700977E-04
      RBC(4,10) = -4.95991087393098E-04     ZBS(4,10) =  2.43875640076056E-05
      RBC(5,10) = -1.94520415280627E-04     ZBS(5,10) =  1.55759001593971E-04
      RBC(6,10) = -6.94143617569942E-05     ZBS(6,10) =  4.40565098025554E-05
    /
    &END

1.  \_\_**Execute the code.**\_\_ \> To simplify execution of the code,
    the VMEC compilation scripts create a directory called \'bin\' in
    your home () directory. Symbolic links are then placed there
    pointing to each of the compiled codes in their respective
    \'Vrelease\' subdirectories. In practice, the screen output from
    VMEC should be redirected to a log file and put in the background
    (\>& log.ncsx\_c09r00\_fixed &) but for this tutorial we will simply
    instruct the user to run the code so they can see how VMEC executes.
    This is done by passing the suffix of the input file to the VMEC
    code through the command line.

<!-- -->

    >~/bin/xvmec2000 ncsx_c09r00_fixed
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SEQ =    1 TIME SLICE  0.0000E+00
      PROCESSING INPUT.ncsx_c09r00_fixed
      THIS IS VMEC2000, A 3D EQUILIBRIUM CODE, VERSION 8.47
      Lambda: Full Radial Mesh. L-Force: hybrid full/half.
      COMPUTER: computer.domain.net   OS: Linux   RELEASE: 2.6.18-194.17.4.el5  DATE = Sep 06,2011  TIME = 12:43:39
      NS =    9 NO. FOURIER MODES =  137 FTOLV =  1.000E-06 NITER =   2500
     ITER    FSQR      FSQZ      FSQL     RAX(v=0)     WMHD
        1  5.24E+01  6.99E+00  1.63E-01  1.608E+00  3.7781E+00
      172  9.30E-07  2.96E-07  2.69E-07  1.610E+00  3.5200E+00
      NS =   28 NO. FOURIER MODES =  137 FTOLV =  1.000E-08 NITER =   2500
     ITER    FSQR      FSQZ      FSQL     RAX(v=0)     WMHD
        1  4.91E-02  2.35E-02  3.37E-04  1.610E+00  3.5199E+00
      200  7.25E-08  2.33E-08  1.02E-08  1.609E+00  3.5196E+00
      300  9.94E-09  3.01E-09  1.35E-09  1.609E+00  3.5196E+00
      NS =   49 NO. FOURIER MODES =  137 FTOLV =  1.000E-11 NITER =   2500
     ITER    FSQR      FSQZ      FSQL     RAX(v=0)     WMHD
        1  1.51E-02  8.86E-03  6.66E-06  1.609E+00  3.5196E+00
      200  5.23E-08  2.49E-08  9.58E-10  1.609E+00  3.5195E+00
      400  2.34E-09  8.33E-10  2.14E-10  1.608E+00  3.5195E+00
      600  3.77E-10  1.11E-10  2.99E-11  1.608E+00  3.5195E+00
      800  7.51E-11  1.87E-11  3.21E-12  1.608E+00  3.5195E+00
      989  9.93E-12  3.22E-12  5.38E-13  1.608E+00  3.5195E+00
     EXECUTION TERMINATED NORMALLY
     FILE : ncsx_c09r00_fixed
     NUMBER OF JACOBIAN RESETS =    4
        TOTAL COMPUTATIONAL TIME              81.34 SECONDS
        TIME TO READ IN DATA                   0.04 SECONDS
        TIME TO WRITE DATA TO WOUT             0.12 SECONDS
        TIME IN EQFORCE                        2.65 SECONDS
        TIME IN FOURIER TRANSFORM             24.27 SECONDS
        TIME IN INVERSE FOURIER XFORM         17.28 SECONDS
        TIME IN FORCES                        17.31 SECONDS
        TIME IN BCOVAR                        12.45 SECONDS
        TIME IN RESIDUE                        0.50 SECONDS
        TIME (REMAINDER) IN FUNCT3D            6.41 SECONDS
                 
    >                                                                   

1.  \_\_**Examine the output.**\_\_ \> For this example four files were
    created
    (<file:jxbout_ncsx_c09r00_fixed.nc>,<file:mercier.ncsx_c09r00_fixed>,<file:threed1.ncsx_c09r00_fixed>,
    and <file:wout_ncsx_c09r00_fixed.nc>). As was mentioned before, each
    file had the suffix of the input file appended to it\'s name. This
    allows multiple runs to be stored in the same directory for
    comparison. The \'jxbout\' file contains values for various
    quantities on a grid throughout the simulation domain. The
    \'mercier\' file contains radial profiles (radial index in VMEC is
    denoted by the variable \'s\') of various quantities. The
    \'threed1\' file can be considered an expanded log file where
    various quantities are calculated which were not output to the
    screen. This file is fairly self explanatory. The \'wout\' file is
    the data file for the run. It contains the Fourier Coefficients for
    the magnetic field along with various quantities. A few packages
    exist to visualize this data and the user is encourage to use these
    as templates for their own visualization routines. \>