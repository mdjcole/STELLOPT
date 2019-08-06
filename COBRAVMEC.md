COBRAVMEC
=========

The COBRAVMEC code calculates ideal ballooning stability for
[VMEC](VMEC) equilibria. \<\<toc\>\>

------------------------------------------------------------------------

### Theory

The Code for Ballooning Rapid
Analysis\<ref\>[R. Sanchez, S. P. Hirshman, J. C. Whitson, and A. S. Ware, \"COBRA: An Optimized Code for Fast Analysis of Ideal Ballooning Stability of Three-Dimensional Magnetic Equilibria.\" J. Comp. Physics 161, 576-588 (2000).](http://www.sciencedirect.com/science/article/pii/S0021999100965148)\</ref\>
(CoBRA) solves for ballooning stability given a VMEC equilibria.
Ballooning instabilites are interchange instabilities which are
localized to regions of poor curvature in toroidal systems. In general
the fourth order ballooning equation can be written [math](math)
\\frac{1}{J}\\partial\_\\theta\\left\[\\frac{\\left(\\nabla
\\beta\\right)\^2}{JB\^\\epsilon}\\partial\_\\theta\\Phi\\right\]+\\left(2p\'\\kappa\_W+\\rho\\omega\^2\\frac{\\left(\\nabla
\\beta\\right)\^2}{B\^2}\\right)\\Phi+2\\kappa\_W\\gamma
p\\nabla\\cdot\\vec{\\xi}=0. [math](math) Neglecting compressional
effects this reduces to [math](math)
\\frac{1}{J}\\frac{\\partial}{\\partial
\\theta}\\left(\\frac{\\left(\\nabla
\\beta\\right)\^2\\partial\\Phi}{B\^2J\\partial\\theta}\\right)+2p\'\\kappa\_W\\Phi+\\rho\\omega\^2\\frac{\\left(\\nabla
\\beta\\right)\^2\\Phi}{B\^2}=0. [math](math) Ideal ballooning stability
can thus be analyzed in terms of a linear second-order differential
equation [math](math)

  --------------------------------------------------- --------
  L\_0\\left(y\\right) + \\lambda R\\left(y\\right)   F = 0,
  --------------------------------------------------- --------

[math](math) where [math](math)
L\_0\\left(y\\right)\\equiv\\frac{d}{dy}\\left\[P\\left(y\\right)\\frac{d}{dy}\\right\]+Q\\left(y\\right)
[math](math) while P(y) and R(y) are positive. Here y is the normalized
length along a field line. Ballooning instability is thus characterized
by negative eigenvalues of this equation (lambda).

------------------------------------------------------------------------

### Compilation

COBRAVMEC is a component of the STELLOPT suite of codes. It is contained
within the \'stellopt.zip\' file. Compilation of the STELLOPT suite is
discussed on the [STELLOPT Compilation Page](STELLOPT Compilation).

------------------------------------------------------------------------

### Input Data Format (pre 4.1)

The COBRAVMEC code takes an in\_cobra text file which controls the
execution of the code. The file has the following format

    5
    5
    test
    0.00  1.58  3.14  4.78  5.46
    0.00  1.58  3.14  4.62  5.00
    5
    20  30  40  50  60

The first line is the number of starting angles in zeta. The second line
is the number of starting angles in theta. The third line is the VMEC
output extension. The fourth line are the starting angles in zeta. The
fifth line contains the starting angles in theta. The sixth line
indicates the number of surfaces on which to calculate ballooning growth
rates. The seventh line indicates the specific surfaces on which to
calculate the ballooning growth rates. It should be noted that in this
example since there are 5 zetas and 5 thetas there will be 25 angles
evaluated.

### Input Data Format (version 4.1)

This version supports non-stellarator symmetric equilibria and the input
format is slightly changed.

    !      INPUT FILE
    !      1st line:   extension of WOUT file
    !      2nd line:   k_w, kth  (No. helical wells; mode number)
    !      3rd line:   l_geom_input, l_tokamak_input
    !      4th line:   nini_zeta (No. initial toroidal angles; or ALPHA labels if L_GEOM_INPUT=F .AND. L_TOKAMAK_INPUT=T)
    !      5th line:   init_zeta_v (vector of initial toroidal angles (or labels), in degrees)
    !      6th line:   nini_theta (No. initial polidal angles; or ALPHA labels if L_GEOM_INPUT=F .AND. L_TOKAMAK_INPUT=F)  
    !      7th line:   init_theta_v (vector of initial poloidal angles (or labels), in degrees)
    !      8th line:   ns_surf (No. surfaces where growth rate is to be computed)
    !      9th line:   ns_v (vector of surfaces where growth rate is to be computed)

------------------------------------------------------------------------

### Execution

The code is executed by passing the in\_cobra input file on the command
line. Screen output may be suppressed by passing \'F\' after the input
file.

    yourmachine:0005> ~/bin/xcobravmec in_cobra.test
    ====================================================

      THIS IS THE (VMEC-BASED) COBRA BALLOONING CODE Version 3.00
     DATE = Sep 27,2011  TIME = 15:32:56

    ----------------------------------------------------------------------------------------------------------------
       NS     FLUX-s     ZT_0        TH_0       GR. RATE    IT   POINTS      XMAX   OK?  SYMM    PRES        MERC
    ----------------------------------------------------------------------------------------------------------------
        2    2.08E-02  0.00E+00    0.00E+00   -5.3585E-02    3     201     1.86E+01  T    0    8.60E-02    1.04E+00
        3    4.17E-02  0.00E+00    0.00E+00   -6.6434E-02    3     201     1.80E+01  T    0    8.58E-02    6.20E-01
        4    6.25E-02  0.00E+00    0.00E+00   -7.7881E-02    3     201     1.75E+01  T    0    8.54E-02    5.28E-01
        5    8.33E-02  0.00E+00    0.00E+00   -8.6155E-02    3     201     1.71E+01  T    0    8.50E-02    4.34E-01
    ----------------------------------------------------------------------------------------------------------------

    ZETA0 =  0.000E+00 THETA0 =  0.000E+00 TIME IN COBRA CODE:  1.02E+00 SEC
    ====================================================

------------------------------------------------------------------------

### Output Data Format

The screen output (if not suppressed) is verbose, containing information
about the runs in tabular format. The code outputs a text file called
cobra\_grate.ext file (where ext is the input extension). For each set
of starting points in theta and zeta it outputs a set of growth rates on
each surface requested. The first line in each set is indicates the
starting position in zeta, theta, and the number of surfaces on which
the growth rate was calculated. Then for each surface the code outputs a
line with the surface index and the growth rate. This is repeated for
each field line followed. And example can be seen below

     0.000E+00 0.000E+00   4
       2 -5.35845793E-02
       3 -6.64339650E-02
       4 -7.78813285E-02
       5 -8.61553492E-02
     3.140E+00 0.000E+00   4
       2 -5.19600197E-02
       3 -6.72031207E-02
       4 -7.57780371E-02
       5 -8.64640971E-02

------------------------------------------------------------------------

### Visualization

For each starting position, the code outputs the growthrate as a
function of flux surface. This can be imported and easily plotted with
many plotting packages.

------------------------------------------------------------------------

### Tutorials

Put links to tutorial pages here.