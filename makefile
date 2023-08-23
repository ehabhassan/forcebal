SHELL = /bin/bash
FC        := gfortran
ld        := $(FC)
moddir    := modules
lib       := iforceballib.a
#FCFLAGS   := -O3 -w -g -fsanitize=address -Wall -fcheck=all -C -ffree-line-length-none
FCFLAGS   := -O3 -w -g -O0 -fbacktrace -Wall -fcheck=all -C -ffree-line-length-none
NETCDFLIB := -L/usr/lib/x86_64-linux-gnu
NETCDFINC := -I/usr/include
NETCDFFLG := -lnetcdff -lnetcdf -ldl -lm
#NETCDFFLG := -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm
ARFLAGS   := -crv
F90OBJS   := spec_kind_mod.f90 forcebal_data_mod.f90 \
	     linear1_mod.f90 spline1_mod.f90 \
	     fluxav_mod.f90 nclass_mod.f90  \
	     read_pro_mod.f0- write_mod.f90 write_netcdf_mod.f90 \
	     forcebal.f90
%.o: %.f90
	$(FC) $(FCFLAGS) $(NETCDFINC) $(NETCDFLIB) $(NETCDFFLG) -c $< -o $@
All: main
main: iforcebal
iforcebal: $(lib)(forcebal.o)
	   ar x $(lib) forcebal.o
	   $(FC) -p -o iforcebal forcebal.o $(lib) $(NETCDFLIB) $(NETCDFINC) $(NETCDFFLG)
	   rm forcebal.o
$(lib)(forcebal.o): $(lib)(spec_kind_mod.o) $(lib)(forcebal_data_mod.o) \
		    $(lib)(linear1_mod.o) $(lib)(spline1_mod.o) \
		    $(lib)(fluxav_mod.o) $(lib)(nclass_mod.o) \
                    $(lib)(read_pro_mod.o) $(lib)(write_mod.o) $(lib)(write_netcdf_mod.o) \
		    $(lib)(forcebal.o)

#$(moddir) :
#	mkdir -p $(moddir)
clean:
	rm -fr $(libobjs) *.o *.mod iforcebal iforceballib.a

