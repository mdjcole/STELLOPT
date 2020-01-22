# LIBSTELL Module
"""

"""
# ===================================================================== #
# ===================================================================== #

#This section is to improve compatability between bython 2.7 and python 3x
from __future__ import absolute_import, division, print_function #, unicode_literals
__metaclass__ = type

import utils as _ut

# ========================================================================== #
# ========================================================================== #

# global string variables used for accessing methods in libstell
s1 = ""
s2 = "mp"
s3 = "_"

# ========================================================================== #
# ========================================================================== #


def setup_libstell(branch=None, path=None, verbose=False):
    """
     Allowing interfacing with this library through difference git branches,
     without messing with environment variables, and allow local user paths
     to pre-compiled libraries
    """
    # ====================================================================== #
    import os as _os
    import sys
    global libstell_branch
    global libstell_path
    global qtCreatorPath
    global use_stellopt
    # ====================================================================== #

    if branch is not None:
        libstell_branch = branch
    if path is not None:
        libstell_path = path
        qtCreatorPath = path
        if verbose and not _os.path.exists(libstell_path):
            print('Error in libstell_path setup:  path does not exist')
        # end if
    # end if

    if 'use_stellopt' not in globals():
        use_stellopt = True
    if 'libstell_branch' not in globals():
        # Default to the release branch of the STELLOPT library
        libstell_branch = 'Release'
    if 'libstell_path' not in globals():
        if use_stellopt:
            if verbose:
                print('Attempting to use the STELLOPT library')
            # end if
            try:
                # Default to the stellopt environment variable
                if 'qtCreatorPath' not in globals():
                    qtCreatorPath = _os.environ["STELLOPT_PATH"]
                # end if
                libstell_path = _os.path.join(qtCreatorPath, 'LIBSTELL', libstell_branch)
            except KeyError:
                if verbose:
                    print('Please set environment variable STELLOPT_PATH or specify the libstell path')
                # end if
                sys.exit()
            # end try
        else:
            # library location in same spot as this module
            libstell_path = _os.path.dirname(__file__)
            if verbose:
                print('Attempting to set a local library: %s'%(libstell_path,))
            # end if

            if 'qtCreatorPath' not in globals():
                qtCreatorPath = libstell_path
            # end if
        # end if
    # end if

    if verbose:
        print('Set libstell library path to: %s'%(libstell_path,))
        print('Set qtCreatorPath for the VMECplot GUI to: %s'%(qtCreatorPath,))
    # end if
# end def


def load_libstell():
    import os as _os
    import sys
    import platform
    import ctypes as ct
#    import numpy.ctypeslib as npct

    global libstell_path
    global qtCreatorPath
    if 'libstell_path' not in globals():
        setup_libstell()
    # end if

    # ====================================================================== #

    os = platform.system()
    is_64bits = sys.maxsize > 2**32

    if is_64bits:   # 64-bit architecture
        # ls_handle = ct.c_longlong
        if os=='Windows':
            libstell =  ct.cdll.LoadLibrary( _os.path.join( libstell_path,"libstell64.dll"))
#            libstell = npct.load_library( _os.path.join( libstell_path,"libstell64.dll"),".")
        elif os=='Linux':
            libstell =  ct.cdll.LoadLibrary( _os.path.join( libstell_path,"libstell64.so"))
#            libstell = npct.load_library( _os.path.join( libstell_path,"libstell64.so"),".")
    else:
        # ls_handle = ct.c_long
        if os=='Windows':
            libstell =  ct.cdll.LoadLibrary( _os.path.join( libstell_path,"libstell.dll"))
#            libstell = npct.load_library( _os.path.join( libstell_path,"libstell.dll"),".")
        elif os=='Linux':
            libstell =  ct.cdll.LoadLibrary( _os.path.join( libstell_path,"libstell.so"))
#            libstell = npct.load_library( _os.path.join( libstell_path,"libstell.so"),".")
        # end if
    # end if

    # ====================================================================== #

    return libstell, libstell_path
# end def load_libstell

# ========================================================================== #


def read_vmec(file, branch=None):
    # Load Libraries
    import ctypes as ct
    import numpy.ctypeslib as npct
#    import numpy as np
    libstell, qtCreatorPath = load_libstell()

    # Read File
    read_wout = getattr(libstell,s1+'read_wout_mod_'+s2+'_readw_and_open'+s3)
    read_wout.argtypes=[ct.c_char_p, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.c_long]
    read_wout.restype=None
    ierr = ct.c_int(0)
    iopen = ct.c_int(0)
    read_wout(file.encode('UTF-8'), ct.byref(ierr), ct.byref(iopen), len(file))
    # Setup Arrays
    vmec_data={}
    # Check
    if not (ierr.value == 0):
        return vmec_data
    # Logical
    varlist=['lasym','lthreed','lwout_opened']
    for temp in varlist:
        vmec_data[temp]=ct.c_bool.in_dll(libstell,s1+'read_wout_mod_'+s2+'_'+temp+s3).value
    # Integers
    varlist=['ns','nfp','mpol','ntor','mnmax','mnmax_nyq','iasym','ierr_vmec']
    for temp in varlist:
        vmec_data[temp]=ct.c_int.in_dll(libstell,s1+'read_wout_mod_'+s2+'_'+temp+s3).value
    # Doubles
    varlist=['wb','wp','gamma','pfac','rmax_surf','rmin_surf','zmax_surf',\
             'aspect','betatot','betapol','betator','betaxis','b0','version_',\
             'ionlarmor','volavgb','fsql','fsqr','fsqz','ftolv','aminor','rmajor',\
             'volume','rbtor','rbtor0','itor','machsq']
    for temp in varlist:
        vmec_data[temp]=ct.c_double.in_dll(libstell,s1+'read_wout_mod_'+s2+'_'+temp+s3).value
    # REAL Arrays (ns)
    varlist = ['iotas','iotaf','presf','phipf','chipf','chi','phi','mass',\
               'pres','beta_vol','phip','buco','bvco','vp','overr','jcuru',\
               'jcurv','specw','jdotb','dmerc','dwell','dcurr','dgeod','equif']
    arr_size = vmec_data['ns']
    ftemp = ct.POINTER(ct.c_double)
    for temp in varlist:
        vmec_data[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'read_wout_mod_'+s2+'_'+temp+s3),(arr_size,1))
    # REAL Arrays (mnmax)
    ftemp = ct.POINTER(ct.c_double)
    vmec_data['xm']=npct.as_array(ftemp.in_dll(libstell,s1+'read_wout_mod_'+s2+'_xm'+s3),(vmec_data['mnmax'],1))
    vmec_data['xn']=npct.as_array(ftemp.in_dll(libstell,s1+'read_wout_mod_'+s2+'_xn'+s3),(vmec_data['mnmax'],1))
    vmec_data['xm_nyq']=npct.as_array(ftemp.in_dll(libstell,s1+'read_wout_mod_'+s2+'_xm_nyq'+s3),(vmec_data['mnmax_nyq'],1))
    vmec_data['xn_nyq']=npct.as_array(ftemp.in_dll(libstell,s1+'read_wout_mod_'+s2+'_xn_nyq'+s3),(vmec_data['mnmax_nyq'],1))

    ## Array values 1D
    ftemp=ct.POINTER(ct.c_double)  # analysis:ignore
    ns = vmec_data['ns']
    mnmax = vmec_data['mnmax']
    mnmax_nyq = vmec_data['mnmax_nyq']

    ## 2D Arrays
    mn2d_size = (ns, mnmax)
    mn2d_nyq_size = (ns, mnmax_nyq)
    fmn=ct.POINTER(ct.c_double)
    for key in ['rmnc','zmns','lmns']:
        vmec_data[key]=npct.as_array(fmn.in_dll(libstell, '%sread_wout_mod_%s_%s%s'%(s1,s2,key,s3)), mn2d_size) #ns,mnmax format
    for key in ['bmnc','gmnc','bsupumnc', 'bsupvmnc', 'bsubsmns', 'bsubumnc', 'bsubvmnc', 'currumnc', 'currvmnc']:
        vmec_data[key]=npct.as_array(fmn.in_dll(libstell, '%sread_wout_mod_%s_%s%s'%(s1,s2,key,s3)), mn2d_nyq_size) #ns,mnmax format
    if vmec_data['iasym']:
        for key in ['rmns','zmnc','lmnc']:
            vmec_data[key]=npct.as_array(fmn.in_dll(libstell, '%sread_wout_mod_%s_%s%s'%(s1,s2,key,s3)), mn2d_size) #ns,mnmax format
        for key in ['bmns','gmns','bsupumns', 'bsupvmns', 'bsubsmnc', 'bsubumns', 'bsubvmns', 'currumns', 'currvmns']:
            vmec_data[key]=npct.as_array(fmn.in_dll(libstell, '%sread_wout_mod_%s_%s%s'%(s1,s2,key,s3)), mn2d_nyq_size) #ns,mnmax format
        # end for
    # end if

#    vmec_data['rmnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_rmnc'+s3),mn2d_size) #ns,mnmax format
#    vmec_data['zmns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_zmns'+s3),mn2d_size) #ns,mnmax format
#    vmec_data['lmns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_lmns'+s3),mn2d_size) #ns,mnmax format
#    vmec_data['bmnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bmnc'+s3),mn2d_nyq_size) #ns,mnmax format
#    vmec_data['gmnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_gmnc'+s3),mn2d_nyq_size) #ns,mnmax format
#    vmec_data['bsupumnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsupumnc'+s3),mn2d_nyq_size) #ns,mnmax format
#    vmec_data['bsupvmnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsupvmnc'+s3),mn2d_nyq_size) #ns,mnmax format
#    vmec_data['bsubsmns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsubsmns'+s3),mn2d_nyq_size) #ns,mnmax format
#    vmec_data['bsubumnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsubumnc'+s3),mn2d_nyq_size) #ns,mnmax format
#    vmec_data['bsubvmnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsubvmnc'+s3),mn2d_nyq_size) #ns,mnmax format
#    vmec_data['currumnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_currumnc'+s3),mn2d_nyq_size) #ns,mnmax format
#    vmec_data['currvmnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_currvmnc'+s3),mn2d_nyq_size) #ns,mnmax format
#    if vmec_data['iasym']:
#        vmec_data['rmns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_rmns'+s3),mn2d_size) #ns,mnmax format
#        vmec_data['zmnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_zmnc'+s3),mn2d_size) #ns,mnmax format
#        vmec_data['lmnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_lmnc'+s3),mn2d_size) #ns,mnmax format
#        vmec_data['bmns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bmns'+s3),mn2d_nyq_size) #ns,mnmax format
#        vmec_data['gmns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_gmns'+s3),mn2d_nyq_size) #ns,mnmax format
#        vmec_data['bsupumns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsupumns'+s3),mn2d_nyq_size) #ns,mnmax format
#        vmec_data['bsupvmns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsupvmns'+s3),mn2d_nyq_size) #ns,mnmax format
#        vmec_data['bsubsmnc']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsubsmnc'+s3),mn2d_nyq_size) #ns,mnmax format
#        vmec_data['bsubumns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsubumns'+s3),mn2d_nyq_size) #ns,mnmax format
#        vmec_data['bsubvmns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_bsubvmns'+s3),mn2d_nyq_size) #ns,mnmax format
#        vmec_data['currumns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_currumns'+s3),mn2d_nyq_size) #ns,mnmax format
#        vmec_data['currvmns']=npct.as_array(fmn.in_dll(libstell,s1+'read_wout_mod_'+s2+'_currvmns'+s3),mn2d_nyq_size) #ns,mnmax format

    # Free memory (don't do this as python accesses this memory)
    #read_wout_dealloc = getattr(libstell,s1+'read_wout_mod_'+s2+'_read_wout_deallocate'+s3)
    #read_wout_dealloc()

    # Correct Arrays (mn-nv) to (mn+nv)
    vmec_data['xn'] = -vmec_data['xn']
    # Put on full grid
    vmec_data['buco'] = h2f(vmec_data['buco'],ns)
    vmec_data['bvco'] = h2f(vmec_data['bvco'],ns)
    vmec_data['vp'] = h2f(vmec_data['vp'],ns)
    vmec_data['overr'] = h2f(vmec_data['overr'],ns)
    vmec_data['specw'] = h2f(vmec_data['specw'],ns)
    # Put matrix quantities on full grid
    for key in ['bmnc','gmnc','lmns','bsupumnc','bsupvmnc','bsubsmns','bsubumnc','bsubvmnc']:
        vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
        vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
        vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]
    if vmec_data['iasym']:
        for key in ['bmns','gmns','lmnc','bsupumns','bsupvmns','bsubsmnc','bsubumns','bsubvmns']:
            vmec_data[key][0,:] = 1.5 * vmec_data[key][1,:] - 0.5 * vmec_data[key][2,:]
            vmec_data[key][1:ns-2,:] = 0.5 * (vmec_data[key][1:ns-2,:] + vmec_data[key][2:ns-1,:])
            vmec_data[key][ns-1,:] = 2.0 * vmec_data[key][ns-2,:] - vmec_data[key][ns-3,:]
    # Return the data dictionary
    return vmec_data

# ========================================================================= #

def h2f(var, ns):
#    import numpy as np
#    temp = np.zeros((ns,1))
#    temp[0] = 1.5 * var[0] - 0.5 * var[1]
#    temp[1:ns-1] = 0.5* (var[0:ns-2] + var[1:ns-1])
#    temp[ns-1] = 1.5 * var[ns-2] - 0.5 * var[ns-3]
#    return temp
    return _ut.h2f(var, ns)

def cfunct(theta, zeta, fmnc, xm, xn):
#    import numpy as np
#    f=0
#    (ns,mn)=fmnc.shape
#    lt = len(theta)
#    lz = len(zeta)
#    mt=np.matmul(xm,theta.T)
#    nz=np.matmul(xn,zeta.T)
#    cosmt=np.cos(mt)
#    sinmt=np.sin(mt)
#    cosnz=np.cos(nz)
#    sinnz=np.sin(nz)
#    f = np.zeros((ns,lt,lz))
#    fmn = np.ndarray((mn,lt))
#    for k in range(ns):
#        fmn = np.broadcast_to(fmnc[k,:],(lt,mn)).T
#        fmncosmt=(fmn*cosmt).T
#        fmnsinmt=(fmn*sinmt).T
#        f[k,:,:]=np.matmul(fmncosmt, cosnz)-np.matmul(fmnsinmt, sinnz)
#    return f
    return _ut.cfunct(theta, zeta, fmnc, xm, xn)

def sfunct(theta, zeta, fmnc, xm, xn):
#    import numpy as np
#    f=0
#    (ns,mn)=fmnc.shape
#    lt = len(theta)
#    lz = len(zeta)
#    mt=np.matmul(xm,theta.T)
#    nz=np.matmul(xn,zeta.T)
#    cosmt=np.cos(mt)
#    sinmt=np.sin(mt)
#    cosnz=np.cos(nz)
#    sinnz=np.sin(nz)
#    f = np.zeros((ns,lt,lz))
#    fmn = np.ndarray((mn,lt))
#    for k in range(ns):
#        fmn = np.broadcast_to(fmnc[k,:],(lt,mn)).T
#        f[k,:,:]=np.matmul((fmn*sinmt).T,cosnz)+np.matmul((fmn*cosmt).T,sinnz)
#    return f
    return _ut.sfunct(theta, zeta, fmnc, xm, xn)

def calc_jll(vmec_data, theta, zeta ):
#    # CALC_JLL(vmec_data,theta,zeta) Calculates the parallel current density.
#    # This funciton takes a VMEC data structure (as read by read_vmec) and
#    # theta/zeta arrays as input and outputs the parallel current density.
#
#    # Example usage (Matlab)
#    #      theta=0:2*pi/359:2*pi;
#    #      zeta=0:2*pi/63:2*pi;
#    #      data=read_vmec('wout.test');        % Reads VMEC wout file
#    #      jll=calc_jll(vmec_data,theta,zeta); % Calculate the current
#    #
#    # Example usage (Python)
#    #      theta=np.linspace(0, 2*np.pi, 360)
#    #      zeta=np.linspace(0, 2*np.pi, 64)
#    #      vmec_data=read_vmec('wout.nc')
#    #      jll=calc_jll(vmec_data, theta, zeta)
#
#
#    # Maintained by: Samuel Lazerson (lazerson@pppl.gov)
#    # Version:       1.00
#
#    b =cfunct(theta,zeta,vmec_data['bmnc'],    vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#    g =cfunct(theta,zeta,vmec_data['gmnc'],    vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#    bu=cfunct(theta,zeta,vmec_data['bsubumnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#    bv=cfunct(theta,zeta,vmec_data['bsubvmnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#    ju=cfunct(theta,zeta,vmec_data['currumnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#    jv=cfunct(theta,zeta,vmec_data['currvmnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#
#    if (vmec_data['iasym']):
#        b =b +sfunct(theta,zeta,vmec_data['bmns'],    vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#        g =g +sfunct(theta,zeta,vmec_data['gmns'],    vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#        bu=bu+sfunct(theta,zeta,vmec_data['bsubumns'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#        bv=bv+sfunct(theta,zeta,vmec_data['bsubvmns'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#        ju=ju+sfunct(theta,zeta,vmec_data['currumns'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#        jv=jv+sfunct(theta,zeta,vmec_data['currvmns'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
#
#
#    jll = (bu*ju+bv*jv)/(g*b)
#    return jll
    return _ut.calc_jll(vmec_data, theta, zeta, use_nyq=True)

# ========================================================================= #
# ========================================================================= #


def torocont(r,z,val,s):
#    import numpy as np
    import matplotlib.pyplot as pyplot
    h=pyplot.axes(xlabel='R [m]',ylabel='Z [m]',aspect='equal')
    pyplot.pcolormesh(r[:,:,s],z[:,:,s],val[:,:,s],cmap='jet',shading='gouraud',axes=h)
    pyplot.show()
    return h

def toroslice(r,zeta,z,s):
    import numpy as np
    import matplotlib.pyplot as pyplot
    h=pyplot.axes(xlabel='R [m]',ylabel='Z [m]',aspect='equal')
    if (s[0] == 0):
        pyplot.plot(r[0,0,zeta],z[0,0,zeta],'+',color='black',axes=h)
        pyplot.plot(np.transpose(r[s[1:],:,zeta]),np.transpose(z[s[1:],:,zeta]),color='black',axes=h)
    else:
        for i in range(0,):
            pyplot.plot(np.transpose(r[s,:,zeta]),np.transpose(z[s,:,zeta]),color='black',axes=h)
    pyplot.show()
    return h

def isotoro(r,z,zeta,svals,*args,**kwargs):
    import numpy as np
    import matplotlib.pyplot as pyplot
#    import mpl_toolkits.mplot3d as mplot3d
    import math as math
    import matplotlib.tri as mtri
    #from mayavi import mlab
    nr = np.size(svals)
    if (nr == 1):
        s= [svals]
        nr = 1
    else:
        s=svals
    nt = np.size(r,1)
    nz = np.size(r,2)
    vertex = np.zeros((nt*nz,3,nr))
    for k in range(0,nr):
        ivertex = 0
        # ifaces = 0
        for j in range(0,nz):
            for i in range(0,nt):
                vertex[ivertex,0,k]=r[s[k],i,j]*math.cos(zeta[j])
                vertex[ivertex,1,k]=r[s[k],i,j]*math.sin(zeta[j])
                vertex[ivertex,2,k]=z[s[k],i,j]
                ivertex = ivertex + 1
    u = np.linspace(0, 1, endpoint=True, num=nt)
    v = np.linspace(0, 1, endpoint=True, num=nz)
    u, v = np.meshgrid(u, v)
    u, v = u.flatten(), v.flatten()
    tri = mtri.Triangulation(u, v)
    test=len(kwargs)
    fig=kwargs.pop('fig',pyplot.figure())
    h=kwargs.pop('axes',fig.add_subplot(111,projection='3d'))
    for k in range(0,nr):
        if (len(args)==0):
            tsurf=h.plot_trisurf(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], triangles=tri.triangles,color='red',shade='yes',linewidths=0.0)
            #tsurf=mlab.triangular_mesh(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], tri.triangless)
        else:
            # Matplotlib way (SLOW)
            vals = args[0][s[k],:,:].T.flatten()
            colors = np.mean(vals[tri.triangles], axis=1)
            tsurf=h.plot_trisurf(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], triangles=tri.triangles,cmap='jet',shade='yes',linewidths=0.0)
            tsurf.set_array(colors)
            tsurf.autoscale()
            #MAYAVI Way (need to figure out how to embed)
            #h    = mlab.figure()
            #vals = args[0][s[k],:,:].T.flatten()
            #tsurf=mlab.triangular_mesh(vertex[:,0,k],vertex[:,1,k],vertex[:,2,k], tri.triangles, scalars=vals, colormap='jet',figure=h)
            #print(type(tsurf))
    if (test==0):
        pyplot.show()
    return h

# ========================================================================= #
# ========================================================================= #

def safe_close(iunit):
    import ctypes as ct
    # Load Libraries
    libstell, qtCreatorPath = load_libstell()

    # Handle interface
    safe_close_h = getattr(libstell,s1+'safe_open_mod_'+s2+'_safe_close'+s3)
    safe_close_h.restype=None
    iunit_temp = ct.c_int(iunit)
    safe_close_h(ct.byref(iunit_temp))
    return

def safe_open(iunit,istat,filename,filestat,fileform,record_in,access_in,delim_in):

    # Load Libraries
    import ctypes as ct
    libstell, qtCreatorPath = load_libstell()

    # Handle interface
    safe_open_h = getattr(libstell,s1+'safe_open_mod_'+s2+'_safe_open'+s3)
    # SUBROUTINE safe_open(int iunit, int istat, char filename, char filestat, char fileform, int record_in, char access_in, char delim_in)
    safe_open_h.argtypes= [ ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.c_char_p, ct.c_char_p, ct.c_char_p, \
        ct.POINTER(ct.c_int), ct.c_char_p, ct.c_char_p, \
        ct.c_long, ct.c_long, ct.c_long, ct.c_long, ct.c_long]
    safe_open_h.restype=None
    iunit_temp = ct.c_int(iunit)
    istat_temp = ct.c_int(istat)
    record_in_temp = ct.c_int(record_in)
    opt1 = ct.c_bool(True)  # analysis:ignore
    opt2 = ct.c_bool(True)  # analysis:ignore
    opt3 = ct.c_bool(True)  # analysis:ignore
    safe_open_h(ct.byref(iunit_temp),ct.byref(istat_temp), \
        filename.encode('UTF-8'),filestat.encode('UTF-8'),fileform.encode('UTF-8'),\
        ct.byref(record_in_temp),access_in.encode('UTF-8'),delim_in.encode('UTF-8'), \
        len(filename),len(filestat),len(fileform),len(access_in),len(delim_in))
    istat = istat_temp
    iunit = iunit_temp
    istat = istat_temp
    return istat

def read_indata_namelist(iunit,istat):
    # Load Libraries
    import ctypes as ct
    import numpy.ctypeslib as npct
    libstell, _ = load_libstell()

    # Handle interface
    read_indata_namelist = getattr(libstell,s1+'vmec_input_'+s2+'_read_indata_namelist'+s3)
    #SUBROUTINE read_indata_namelist (iunit, istat)
    read_indata_namelist.argtypes = [ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    read_indata_namelist.restype=None
    iunit_temp = ct.c_int(iunit)
    istat_temp = ct.c_int(istat)
    read_indata_namelist(ct.byref(iunit_temp),ct.byref(istat_temp))
    istat = istat_temp
    iunit = iunit_temp
    # Setup Arrays
    indata_namelist={}
    # Logicals
    varlist=['lpofr','lmac','lfreeb','lrecon','loldout','ledge_dump','lasym','lforbal','lrfp',\
             'lmovie','lmove_axis','lwouttxt','ldiagno','lmoreiter','lfull3d1out','l_v3fit',\
             'lspectrum_dump','loptim','lgiveup','lbsubs','lgiveup']
    for temp in varlist:
        indata_namelist[temp]=ct.c_bool.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3).value
    # Integers
    varlist=['nfp','ncurr','nsin','niter','nstep','nvacskip','mpol','ntor','ntheta','nzeta', \
             'mfilter_fbdy','nfilter_fbdy','max_main_iterations','omp_num_threads',\
             'imse','isnodes','itse','ipnodes','iopt_raxis','imatch_phiedge','nflxs']
    for temp in varlist:
        indata_namelist[temp]=ct.c_int.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3).value
    # Reals
    varlist=['time_slice','curtor','delt','ftol','tcon0','gamma','phiedge','phidiam',\
             'sigma_current','sigma_delphid','tensi','tensp','tensi2','fpolyi','presfac',\
             'mseangle_offset','pres_offset','mseangle_offsetm','spres_ped','bloat',\
             'pres_scale','prec2d_threshold','bcrit','fgiveup']
    for temp in varlist:
        indata_namelist[temp]=ct.c_double.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3).value
    # Get integers defined elsewhere (Hardcode for now, not sure how to get them)
    #indata_namelist['nbsetsp']=ct.c_int.in_dll(libstell,s1+'vsvd0_'+s2+'_nbsetsp').value
    # Integer Arrays (100)
    varlist = ['ns_array','niter_array']
    arr_size=100
    ftemp = ct.c_int*arr_size
    for temp in varlist:
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
    # Note that we skip some arrays related to recon stuff because we don't need them and we
    # need to figure out how to pull stuff from other modules see the above issue.
    # Real 2D Arrays (ntord=101,mpol1d=100)
    varlist = ['rbs','zbc','rbc','zbs']
    arr_size1=2*101+1
    arr_size2=100+1
    ftemp = ct.c_double*arr_size1*arr_size2
    for temp in varlist:
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size1,arr_size2))
    # REAL Arrays (21)
    varlist = ['am','ai','ac','ah','at']
    arr_size=21
    ftemp = ct.c_double*arr_size
    for temp in varlist:
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
    # REAL Arrays (20)
    varlist = ['aphi']
    arr_size=20
    ftemp = ct.c_double*arr_size
    for temp in varlist:
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
    # REAL Arrays (ndatafmax=101)
    varlist = ['am_aux_s','am_aux_f','ac_aux_s','ac_aux_f','ai_aux_s','ai_aux_f',\
               'ah_aux_s','ah_aux_f','at_aux_s','at_aux_f']
    arr_size=101
    ftemp = ct.c_double*arr_size
    for temp in varlist:
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
    # REAL Arrays (ntord+1=102)
    varlist = ['raxis','zaxis','raxis_cc','raxis_cs','zaxis_cc','zaxis_cs']
    arr_size=102
    ftemp = ct.c_double*arr_size
    for temp in varlist:
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
    # REAL Arrays (100)
    varlist = ['ftol_array']
    arr_size=100
    ftemp = ct.c_double*arr_size
    for temp in varlist:
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
    # REAL Arrays (nigroup=300)
    varlist = ['extcur']
    arr_size=300
    ftemp = ct.c_double*arr_size
    for temp in varlist:
        indata_namelist[temp]=npct.as_array(ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3),(arr_size,1))
    # Charater arrays
    varlist = ['pcurr_type','piota_type','pmass_type','pt_type','ph_type']
    arr_size=20
    ftemp = ct.c_char*arr_size
    for temp in varlist:
        indata_namelist[temp]=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_'+temp+s3).value.decode('UTF-8')
    ftemp = ct.c_char*200
    indata_namelist['mgrid_file']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_mgrid_file'+s3).value.decode('UTF-8')
    ftemp = ct.c_char*200
    indata_namelist['trip3d_file']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_trip3d_file'+s3).value.decode('UTF-8')
    ftemp = ct.c_char*10
    indata_namelist['precon_type']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_precon_type'+s3).value.decode('UTF-8')
    ftemp = ct.c_char*120
    indata_namelist['arg1']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_arg1'+s3).value.decode('UTF-8')
    ftemp = ct.c_char*100
    indata_namelist['input_extension']=ftemp.in_dll(libstell,s1+'vmec_input_'+s2+'_input_extension'+s3).value.decode('UTF-8')
    return indata_namelist

def set_module_var(module,var,val):
    # Load Libraries
    import numpy as np
    import ctypes as ct
    libstell, qtCreatorPath = load_libstell()

    if type(val) == bool:
        f = ct.c_bool
    elif type(val) == int:
        f = ct.c_int
    elif type(val) == float:
        f = ct.c_double
    elif type(val) == str:
        n = len(val)
        f = ct.c_char*n
    elif type(val) == np.ndarray:
        if type(val[0]) == bool:
            tt = ct.c_bool
        elif type(val[0]) == np.int32:
            tt = ct.c_int
        elif type(val[0]) == np.float64:
            tt = ct.c_double
        else:
            print('   Unrecognized type:',type(val[0]))
            return
        n = val.ndim
        f = tt*val.size
        #print(n,val.size)
    else:
        print('   Unrecognized type:',type(val))
        return
    temp=f.in_dll(libstell,s1+''+module+'_'+s2+'_'+var)
    if type(val) == np.ndarray:
        if n==1:
            for i,col in enumerate(val):
                temp[i] = val[i]
    elif type(val) == str:
        temp.value = val.encode('UTF-8')
    else:
        temp.value = val
    return

def write_indata_namelist(iunit,istat):
    # Load Libraries
    import ctypes as ct
    libstell, qtCreatorPath = load_libstell()

    # Handle interface
    write_indata_namelist = getattr(libstell,s1+'vmec_input_'+s2+'_write_indata_namelist'+s3)
    #SUBROUTINE read_indata_namelist (iunit, istat)
    write_indata_namelist.argtypes = [ct.POINTER(ct.c_int),ct.POINTER(ct.c_int)]
    write_indata_namelist.restype=None
    iunit_temp = ct.c_int(iunit)
    istat_temp = ct.c_int(istat)
    write_indata_namelist(ct.byref(iunit_temp),ct.byref(istat_temp))
    #istat = istat_temp
    #iunit = iunit_temp
    return

def pcurr(xx):
    # Load Libraries
    import ctypes as ct
    libstell, qtCreatorPath = load_libstell()

    pcurr_func = getattr(libstell,'pcurr_')
    #SUBROUTINE pcurr (xx)
    pcurr_func.argtypes = [ct.POINTER(ct.c_double)]
    pcurr_func.restype=ct.c_double
    xx_temp = ct.c_double(xx)
    val = pcurr_func(ct.byref(xx_temp))
    return val;

def pmass(xx):
    # Load Libraries
    import ctypes as ct
    libstell, qtCreatorPath = load_libstell()

    pmass_func = getattr(libstell,'pmass_')
    #SUBROUTINE piota (xx)
    pmass_func.argtypes = [ct.POINTER(ct.c_double)]
    pmass_func.restype=ct.c_double
    xx_temp = ct.c_double(xx)
    val = pmass_func(ct.byref(xx_temp))
    return val;

def piota(xx):
    # Load Libraries
    import ctypes as ct
    libstell, qtCreatorPath = load_libstell()

    piota_func = getattr(libstell,'piota_')
    #SUBROUTINE piota (xx)
    piota_func.argtypes = [ct.POINTER(ct.c_double)]
    piota_func.restype=ct.c_double
    xx_temp = ct.c_double(xx)
    val = piota_func(ct.byref(xx_temp))
    return val;



# ========================================================================== #
# ========================================================================== #








