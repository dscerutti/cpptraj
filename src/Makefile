# Cpptraj standalone Makefile
include ../config.h

include cpptrajfiles

OBJECTS=$(SOURCES:.cpp=.o) $(CSOURCES:.c=.o)

include ambpdbfiles

AMBPDB_OBJECTS=$(AMBPDBSOURCES:.cpp=.o)

CUDA_SOURCES=cuda_kernels/core_kernels.cu cuda_kernels/kernel_wrappers.cu

# Default target: cpptraj only
all: cpptraj$(SFX)$(EXE)

# Install all targets
install: $(INSTALL_TARGETS) 

# Cpptraj ------------------------------
install_cpptraj: cpptraj$(SFX)$(EXE)
	mv cpptraj$(SFX)$(EXE) $(CPPTRAJBIN)

cpptraj$(SFX)$(EXE): $(OBJECTS) $(FFT_TARGET) $(READLINE_TARGET) $(CUDA_TARGET) $(XDRFILE_TARGET) $(ARPACK_TARGET)
	$(CXX) -o cpptraj$(SFX)$(EXE) $(OBJECTS) $(CUDA_TARGET) $(FFT_TARGET) $(READLINE_LIB) $(CPPTRAJ_LIB) $(LDFLAGS)

# AmbPDB -------------------------------
install_ambpdb: ambpdb$(EXE)
	mv ambpdb$(EXE) $(CPPTRAJBIN)

ambpdb$(EXE): $(AMBPDB_OBJECTS) $(XDRFILE_TARGET)
	$(CXX) -o ambpdb$(EXE) $(AMBPDB_OBJECTS) $(LDFLAGS)

# libcpptraj ---------------------------
# Rule to make libcpptraj-specific objects
%.LIBCPPTRAJ.o : %.cpp
	$(CXX) $(DIRECTIVES) -DLIBCPPTRAJ $(INCLUDE) $(CXXFLAGS) -c -o $@ $<

libcpptraj: $(LIBCPPTRAJ_TARGET)

$(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX): $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(XDRFILE_TARGET) $(ARPACK_TARGET)
	$(CXX) -shared -o $(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX) $(LIBCPPTRAJ_OBJECTS) $(FFT_TARGET) $(CPPTRAJ_LIB) $(LDFLAGS)

nolibcpptraj:
	@echo ""
	@echo "Error: Cannot build libcpptraj; re-configure with '-shared'"
	@echo ""
	@exit 1

# nproc --------------------------------
nproc:
	cd ../util && $(MAKE)

# Internal object/library targets
pub_fft.o:  pub_fft.F90
	$(FC) $(FPPFLAGS) -c $(FFLAGS) -o $@ pub_fft.F90

readline/libreadline.a: readline/Makefile
	cd readline && $(MAKE) all

xdrfile/libxdrfile.a: xdrfile/Makefile
	cd xdrfile && $(MAKE) all

arpack/libarpack.a: arpack/Makefile
	cd arpack && $(MAKE) install

noreadline:
	@echo ""
	@echo "Skipping bundled READLINE build"
	@echo ""

noxdrfile:
	@echo ""
	@echo "Skipping bundled XDRFILE build"
	@echo ""

noarpack:
	@echo ""
	@echo "Skipping bundled ARPACK build"
	@echo ""

cuda_kernels/libcpptraj_cuda.a: $(CUDA_SOURCES)
	cd cuda_kernels && $(MAKE) all

# Dependency targets
findDepend: FindDepend.o 
	$(CXX) -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCES) $(CSOURCES) AmbPDB.cpp > cpptrajdepend

dependclean:
	/bin/rm -f FindDepend.o findDepend

# Clean/uninstall targets
clean:
	/bin/rm -f $(OBJECTS) $(FFT_TARGET) cpptraj$(SFX) AmbPDB.o ambpdb *.LIBCPPTRAJ.o libcpptraj$(SHARED_SUFFIX)
	cd readline && $(MAKE) clean
	cd xdrfile && $(MAKE) clean
	cd arpack && $(MAKE) clean
	cd cuda_kernels && $(MAKE) clean

uninstall:
	/bin/rm -f $(CPPTRAJBIN)/cpptraj$(SFX)$(EXE)
	/bin/rm -f $(CPPTRAJBIN)/ambpdb$(EXE)
	/bin/rm -f $(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX)
	cd readline && make uninstall
	cd xdrfile && make uninstall
	cd cuda_kernels && make uninstall

# Header dependencies
include cpptrajdepend
