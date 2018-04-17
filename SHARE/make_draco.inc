#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(HOME)/bin
  PRECOMP:= fpp -C -P -DLINUX
  COMPILE = ifort
  COMPILE_FREE = ifort
  LINK    = ftn -Bstatic -o
  LINK_AR = ar -ruv
  LINK_C  = $(CC) -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -I$(MKL_HOME)/include -O2 -fp-model strict -xAVX -ip \
            -assume noold_unit_star 
  FLAGS_D = -g -I$(MKL_HOME)/include
  LIBS    = -Wl,-rpath,$(MKL_HOME)/lib/intel64 \
            -L$(MKL_HOME)/lib/intel64 -lmkl_scalapack_lp64 \
            -lmkl_intel_lp64 -lmkl_core -lmkl_sequential \
            -lmkl_blacs_intelmpi_lp64  -lpthread -lm -lpthread

#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpiifort
  MPI_COMPILE_FREE = mpiifort
  MPI_COMPILE_C = 
  MPI_LINK = mpiifort

#######################################################################
#            NAG Options
#######################################################################
  LNAG = T
  NAG_LIB = $(NAGFLIB) $(NAGMKLLIB)

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = -I$(NETCDF_HOME)/include 
  NETCDF_LIB = -Wl,-rpath,$(NETCDF_HOME)/lib -L$(NETCDF_HOME)/lib \
               -lnetcdf -lnetcdff 

#######################################################################
#            NTCC Options
#######################################################################
  LNTCC = F
  NTCC_INC = -I$(NTCCHOME)/mod
  NTCC_LIB = -L$(NTCCHOME)/lib -laladdinsub -ladpak -lcppsub \
             -lcomput -lportlib 

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  HDF5_INC = -I$(HDF5_HOME)/include 
  HDF5_LIB = -Wl,-rpath,$(HDF5_HOME)/lib -L$(HDF5_HOME)/lib \
             -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran 

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC =
  PGPLOT_LIB = 

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = -I$(SILO_DIR)/include
  SILO_LIB = -L$(SILO_DIR)/lib -lsilo

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = T
  LNEO  = T

#######################################################################
#            GENE Options
#######################################################################
  LGENE = F
  GENE_INC = -I$(GENE_PATH)
  GENE_DIR = $(GENE_PATH)
  LIB_GENE = libgene.a
  GENE_LIB = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f -Xlinker \
             -rpath=$(FFTW_HOME)/lib \
             -Wl,-rpath,$(SLEPC_HOME)/intel14.0-mpi.ibm1308-double-cplx/lib \
             -L$(SLEPC_HOME)/intel14.0-mpi.ibm1308-double-cplx/lib -lslepc \
             -Wl,-rpath,$(PETSC_HOME)/lib -L$(PETSC_HOME)/lib -lpetsc \
             -Wl,-rpath,$(MKL_HOME)/lib/intel64 -L$(MKL_HOME)/lib/intel64 \
             -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core \
             -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm 

#######################################################################
#            COILOPT++ Options
#######################################################################
  LCOILOPT = F
  COILOPT_INC = -I$(COILOPT_PATH)
  COILOPTPP_DIR = $(COILOPT_PATH)
  LIB_COILOPTPP = libcoilopt++.a
  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) 

#######################################################################
#            TERPSICHORE Options
#######################################################################
  LTERPSICHORE= F
  TERPSICHORE_INC = -I$(TERPSICHORE_PATH)
  TERPSICHORE_DIR = $(TERPSICHORE_PATH)
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)

#######################################################################
#            TRAVIS Options
#######################################################################
  LTRAVIS= T
  TRAVIS_DIR = $(TRAVIS_PATH)
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB = $(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
               $(TRAVIS_DIR)/mag_conf/lib/$(LIB_MCONF) -lstdc++

#######################################################################
#            REGCOIL Options
#######################################################################
  LREGCOIL= F
  REGCOIL_DIR = $(REGCOIL_PATH)
  REGCOIL_INC = -I$(REGCOIL_DIR) 
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
LIB_SHARE = 

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#   DO NOT EDIT BELOW HERE DO NOT EDIT BELOW HERE

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

# This is where we find the cawk script
#HOME_BIN = $(MYHOME)
HOME_BIN = ../../
MOD1_PATH = -I$(MYHOME)/libstell_dir
MOD_PATH = 
CLEAN_ARCH = 
LIB_DIR = ../../LIBSTELL
LIB = libstell.a

# Handle Which Code is being CALLED
ifeq ($(MAKECMDGOALS),xstelloptv2)
  PRECOMP += -DBNORM -DBEAMS3D_OPT -DJINV_OPT -DTXPORT_OPT -DDIAGNO_OPT -DSKS  -DSKS2
  VMEC_DIR ?= ../../VMEC2000
  LIB_VMEC= libvmec.a
  MOD_PATH+= -I$(VMEC_DIR)/$(LOCTYPE)
  BEAMS3D_DIR = ../../BEAMS3D
  LIB_BEAMS3D = libbeams3d.a
  MOD_PATH+= -I$(BEAMS3D_DIR)/$(LOCTYPE)
  BOOTSJ_DIR = ../../BOOTSJ
  LIB_BOOTSJ = libbootsj.a
  MOD_PATH+= -I$(BOOTSJ_DIR)/$(LOCTYPE)
  BNORM_DIR = ../../BNORM
  LIB_BNORM = libbnorm.a
  MOD_PATH+= -I$(BNORM_DIR)/$(LOCTYPE)
  BOOZ_DIR = ../../BOOZ_XFORM
  LIB_BOOZ= libbooz.a
  MOD_PATH+= -I$(BOOZ_DIR)/$(LOCTYPE)
  COBRA_DIR = ../../COBRAVMEC
  LIB_COBRA = libcobra.a
  MOD_PATH+= -I$(COBRA_DIR)/$(LOCTYPE)
  DIAGNO_DIR = ../../DIAGNO
  LIB_DIAGNO= libdiagno.a
  MOD_PATH+= -I$(DIAGNO_DIR)/$(LOCTYPE)
  JINV_DIR = ../../J_INVARIANT
  LIB_JINV = libjinv.a
  MOD_PATH+= -I$(JINV_DIR)/$(LOCTYPE)
  MGRID_DIR= ../../MAKEGRID
  LIB_MGRID= libmakegrid.a
  MOD_PATH+= -I$(MGRID_DIR)/$(LOCTYPE)
  LIB_LINK=  $(VMEC_DIR)/$(LOCTYPE)/$(LIB_VMEC)  \
	$(BEAMS3D_DIR)/$(LOCTYPE)/$(LIB_BEAMS3D) \
	$(BOOTSJ_DIR)/$(LOCTYPE)/$(LIB_BOOTSJ) \
	$(BNORM_DIR)/$(LOCTYPE)/$(LIB_BNORM) \
	$(BOOZ_DIR)/$(LOCTYPE)/$(LIB_BOOZ) \
	$(COBRA_DIR)/$(LOCTYPE)/$(LIB_COBRA) \
	$(DIAGNO_DIR)/$(LOCTYPE)/$(LIB_DIAGNO) \
	$(JINV_DIR)/$(LOCTYPE)/$(LIB_JINV) \
	$(MGRID_DIR)/$(LOCTYPE)/$(LIB_MGRID) 
  ifeq ($(LDKES),T)
    DKES_DIR = ../../DKES
    LIB_DKES = libdkes.a
    MOD_PATH += -I$(DKES_DIR)/$(LOCTYPE)
    PRECOMP += -DDKES_OPT
    LIB_LINK += $(DKES_DIR)/$(LOCTYPE)/$(LIB_DKES)
  else
    DKES_DIR = 
    LIB_DKES = 
  endif
  ifeq ($(LNEO),T)
    NEO_DIR = ../../NEO
    LIB_NEO = libneo.a
    MOD_PATH += -I$(NEO_DIR)/$(LOCTYPE)
    PRECOMP += -DNEO_OPT
    LIB_LINK += $(NEO_DIR)/$(LOCTYPE)/$(LIB_NEO)
  else
    NEO_DIR = 
    LIB_NEO = 
  endif
else
  LGENE = F
  LCOILOPT = F
  LTERPSICHORE = F
  LTRAVIS = F
endif
ifeq ($(MAKECMDGOALS),xanimec)
  PRECOMP += -D_ANIMEC
endif
ifeq ($(MAKECMDGOALS),xsiesta)
  PRECOMP += -DSKS  -DSKS2
endif
ifeq ($(MAKECMDGOALS),xgtovmi)
  ifeq ($(LPGPLOT),T)
    MOD1_PATH += $(PGPLOT_INC)
    LIB_LINK += $(PGPLOT_LIB)
  endif
endif
ifeq ($(MAKECMDGOALS),xdes_plot)
  ifeq ($(LPGPLOT),T)
    MOD1_PATH += $(PGPLOT_INC)
    LIB_LINK += $(PGPLOT_LIB)
  endif
endif
ifeq ($(MAKECMDGOALS),xvmec2000)
  PRECOMP += -DSKS  -DSKS2
endif

# Other Codes
ifeq ($(LGENE),T)
  MOD1_PATH += $(GENE_INC)
  PRECOMP += -DGENE
  MOD_PATH += -I$(GENE_DIR)
else
  GENE_LIB = 
  GENE_DIR = 
  LIB_GENE = 
endif
ifeq ($(LCOILOPT),T)
  MOD1_PATH += $(COILOPT_INC)
  PRECOMP += -DCOILOPTPP
  MOD_PATH += -I$(COILOPTPP_DIR)
else
  COILOPT_LIB = 
  COILOPTPP_DIR = 
  LIB_COILOPTPP = 
endif
ifeq ($(LTERPSICHORE),T)
  MOD1_PATH += $(TERPSICHORE_INC)
  PRECOMP += -DTERPSICHORE
  MOD_PATH += -I$(TERPSICHORE_DIR)
else
  TERPSICHORE_LIB = 
  TERPSICHORE_DIR = 
  LIB_TERPSICHORE = 
endif
ifeq ($(LTRAVIS),T)
  PRECOMP += -DTRAVIS
else
  TRAVIS_LIB = 
  TRAVIS_DIR = 
  LIB_TRAVIS = 
endif
ifeq ($(LREGCOIL),T)
  MOD1_PATH += $(REGCOIL_INC)
  PRECOMP += -DREGCOIL
  MOD_PATH += -I$(REGCOIL_DIR)
else
  REGCOIL_LIB = 
  REGCOIL_DIR = 
  LIB_REGCOIL = 
endif

# Handle other codes
LIB_LINK += $(GENE_LIB) $(COILOPT_LIB) $(TERPSICHORE_LIB) $(TRAVIS_LIB) $(REGCOIL_LIB)

# Add libraries to LIB_LINK
#LIB_LINK+= $(HOME_BIN)/libstell.a $(LIBS)
LIB_LINK+= $(LIB_DIR)/$(LOCTYPE)/$(LIB) $(LIBS)

# 3rd Party Libs
ifeq ($(LMPI),T)
  PRECOMP += -DMPI_OPT
  VMEC_DIR = ../../VMEC2000
  COMPILE = $(MPI_COMPILE)
  COMPILE_FREE = $(MPI_COMPILE_FREE)
  LINK    = $(MPI_COMPILE) $(FLAGS) -o 
  LINK_C  = $(MPI_COMPILE_C) -shared -Wl,-z,defs
endif
ifeq ($(LNETCDF),T)
  MOD1_PATH += $(NETCDF_INC)
  PRECOMP += -DNETCDF
  LIB_LINK += $(NETCDF_LIB)
  MOD_PATH += $(NETCDF_INC)
endif
ifeq ($(LHDF5),T)
  MOD1_PATH += $(HDF5_INC)
  PRECOMP += -DLHDF5 -DH5_USE_16_API -DHDF5_PAR
  LIB_LINK += $(HDF5_LIB)
  MOD_PATH += $(HDF5_INC)
endif
ifeq ($(LNAG),T)
  PRECOMP += -DNAG
  LIB_LINK += $(NAG_LIB)
endif
ifeq ($(LNTCC),T)
  MOD1_PATH += $(NTCC_INC)
  PRECOMP += -DLADAS -DNTCC
  LIB_LINK += $(NTCC_LIB)
  MOD_PATH += $(NTCC_INC)
endif
ifeq ($(LSILO),T)
  MOD1_PATH += $(SILO_INC)
  PRECOMP += -DSILO_AVAI
  LIB_LINK += $(SILO_LIB)
endif

# This is hear because of MPI above resetting LINK
ifeq ($(MAKECMDGOALS),libstell.a)
  LINK = $(LINK_AR)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif
ifeq ($(MAKECMDGOALS),libstell.so)
  LIB_LINK = $(LIB_SHARE)
  COMPILE += -fPIC
  COMPILE_FREE += -fPIC
  INC_PATH = 
endif

#This must come here
MOD_PATH+= -I