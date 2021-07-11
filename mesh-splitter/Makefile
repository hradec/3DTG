

DOCKER_ARGS?=

PIPE_LIBS_ROOT=/atomo/pipeline/libs/linux/x86_64/pipevfx.5.0.0
CORTEX_VERSION=10.0.0.a82
BOOST_VERSION=1.61.0
GCC_VERSION=$(shell grep gcc ${PIPE_LIBS_ROOT}/cortex/${CORTEX_VERSION}/.build.noBaseLib-boost.${BOOST_VERSION}.done.depend | awk -F'/' '{print $$(NF)}' | sort -V | tail -1)
ILMBASE_VERSION=$(shell grep ilmbase ${PIPE_LIBS_ROOT}/cortex/${CORTEX_VERSION}/.build.noBaseLib-boost.${BOOST_VERSION}.done.depend | awk -F'/' '{print $$(NF)}' | sort -V | tail -1)
OPENEXR_VERSION=$(shell grep openexr ${PIPE_LIBS_ROOT}/cortex/${CORTEX_VERSION}/.build.noBaseLib-boost.${BOOST_VERSION}.done.depend | awk -F'/' '{print $$(NF)}' | sort -V | tail -1)
TBB_VERSION=$(shell grep tbb ${PIPE_LIBS_ROOT}/cortex/${CORTEX_VERSION}/.build.noBaseLib-boost.${BOOST_VERSION}.done.depend | awk -F'/' '{print $$(NF)}' | sort -V | tail -1)
BLOSC_VERSION=$(shell grep blosc ${PIPE_LIBS_ROOT}/cortex/${CORTEX_VERSION}/.build.noBaseLib-boost.${BOOST_VERSION}.done.depend | awk -F'/' '{print $$(NF)}' | sort -V | tail -1)
FREETYPE_VERSION=$(shell grep freetype ${PIPE_LIBS_ROOT}/cortex/${CORTEX_VERSION}/.build.noBaseLib-boost.${BOOST_VERSION}.done.depend | awk -F'/' '{print $$(NF)}' | sort -V | tail -1)
ALEMBIC_VERSION=$(shell grep alembic  ${PIPE_LIBS_ROOT}/cortex/${CORTEX_VERSION}/.build.noBaseLib-boost.${BOOST_VERSION}.done.depend | awk -F'/' '{print $$(NF)}' | sort -V | tail -1)
MPFR_VERSION=$(shell grep mpfr ${PIPE_LIBS_ROOT}/gcc/${GCC_VERSION}/.build.noBaseLib.done.depend | awk -F'/' '{print $$(NF)}' | sort -V | tail -1)
PYTHON_VERSION_MAJOR=$(shell grep python ${PIPE_LIBS_ROOT}/cortex/${CORTEX_VERSION}/.build.noBaseLib-boost.${BOOST_VERSION}.done.depend | awk -F'/' '{print $$(NF)}' | sort -V | tail -1)

ALEMBIC_VERSION=1.7.1

# gcc
ROOT_CXX  		:= ${PIPE_LIBS_ROOT}/gcc/${GCC_VERSION}
CXX 	  		:= ${ROOT_CXX}/bin/g++
CXXFLAGS  		:= -Wno-deprecated -pipe -DBOOST_DISABLE_ASSERTS -fpermissive
CXXFLAGS		+= -ggdb -DDEBUG
# CXXFLAGS		+= -DNDEBUG
CXXFLAGS		+= -DAtUInt=AtUInt32 -D_GLIBCXX_USE_CXX11_ABI=0 -fPIC -I. -nostdinc -nostdinc++
CXXFLAGS		+= -I${ROOT_CXX}/include/
CXXFLAGS		+= -I${ROOT_CXX}/include/c++/
CXXFLAGS		+= -I${ROOT_CXX}/include/c++/${GCC_VERSION}/
# CXXFLAGS		+= -I${ROOT_CXX}/include/c++/${GCC_VERSION}/parallel
CXXFLAGS		+= -I${ROOT_CXX}/include/c++/${GCC_VERSION}/x86_64-pc-linux-gnu
CXXFLAGS		+= -I${ROOT_CXX}/include/c++/${GCC_VERSION}/x86_64-pc-linux-gnu/bits/
CXXFLAGS		+= -I${ROOT_CXX}/lib/gcc/x86_64-pc-linux-gnu/6.3.1/include/
CXXFLAGS		+= -I${ROOT_CXX}/lib/gcc/x86_64-pc-linux-gnu/6.3.1/include/ssp/
CXXFLAGS		+= -I${ROOT_CXX}/lib/gcc/x86_64-pc-linux-gnu/6.3.1/include-fixed/
CXXFLAGS		+= -I/usr/include/
INSTALL_LIBS	+= $(patsubst %, ${ROOT_CXX}/lib/%*.so*, libstdc++ libgcc_s)
LD_LIBRARY_PATH := ${ROOT_CXX}/lib64/:${LD_LIBRARY_PATH}

# mpfr (needed by gcc)
ROOT_MPFR 		:= ${PIPE_LIBS_ROOT}/mpfr/${MPFR_VERSION}
LIBS_MPFR 		+= mpfr
CXXFLAGS 		+= -I${ROOT_MPFR}/include
LD_LIBRARY_PATH := ${ROOT_MPFR}/lib:${LD_LIBRARY_PATH}
LIBS	 		+= $(patsubst %, -l%, $(LIBS_MPFR))
INSTALL_LIBS 	+= $(patsubst %, ${ROOT_MPFR}/lib/lib%*.so*, $(LIBS_MPFR))

# cortexVFX
ROOT_CORTEX		:= ${PIPE_LIBS_ROOT}/cortex/${CORTEX_VERSION}
LIBS_CORTEX		:= IECoreImage IECoreScene IECore IECoreAlembic
CXXFLAGS 		+= -I${ROOT_CORTEX}/include
LIBS	 		+= $(patsubst %, -l%, $(LIBS_CORTEX))
INSTALL_LIBS 	+= $(patsubst %, ${ROOT_CORTEX}/lib/boost${BOOST_VERSION}/lib%*.so, $(LIBS_CORTEX))

# boost
ROOT_BOOST		:= ${PIPE_LIBS_ROOT}/boost/${BOOST_VERSION}
LIBS_BOOST		:= boost_system
CXXFLAGS 		+= -I${ROOT_BOOST}/include
LDFLAGS  		+= -L${ROOT_BOOST}/lib/python${PYTHON_VERSION_MAJOR}
LIBS	 		+= $(patsubst %, -l%, $(LIBS_BOOST))
INSTALL_LIBS 	+= $(wildcard ${ROOT_BOOST}/lib/lib*.so*)

# blosc
ROOT_BLOSC		:= ${PIPE_LIBS_ROOT}/blosc/${BLOSC_VERSION}
CXXFLAGS 		+= -I${ROOT_BLOSC}/include
LDFLAGS  		+= -L${ROOT_BLOSC}/lib/
INSTALL_LIBS 	+= $(wildcard ${ROOT_BLOSC}/lib/lib*.so*)

# freetype
ROOT_FREETYPE	:= ${PIPE_LIBS_ROOT}/freetype/${FREETYPE_VERSION}
CXXFLAGS 		+= -I${ROOT_FREETYPE}/include
LDFLAGS  		+= -L${ROOT_FREETYPE}/lib/
INSTALL_LIBS 	+= $(wildcard ${ROOT_FREETYPE}/lib/lib*.so*)

# openexr
ROOT_OPENEXR	:= ${PIPE_LIBS_ROOT}/openexr/${OPENEXR_VERSION}/boost.${BOOST_VERSION}
LIBS_OPENEXER	+= IlmImf
CXXFLAGS 		+= -I${ROOT_OPENEXR}/include
CXXFLAGS 		+= -I${ROOT_OPENEXR}/include/OpenEXR
LIBS	 		+= $(patsubst %, -l%, $(LIBS_OPENEXER))
INSTALL_LIBS 	+= $(patsubst %, ${ROOT_OPENEXR}/lib/lib%*.so*, $(LIBS_OPENEXER))

# ilmbase
ROOT_ILMBASE	:= ${PIPE_LIBS_ROOT}/ilmbase/${ILMBASE_VERSION}/boost.${BOOST_VERSION}/
LIBS_ILMBASE	:= Half Iex IexMath IlmThread Imath
CXXFLAGS 		+= -I${ROOT_ILMBASE}/include
CXXFLAGS 		+= -I${ROOT_ILMBASE}/include/OpenEXR
LIBS	 		+= $(patsubst %, -l%, $(LIBS_ILMBASE))
INSTALL_LIBS 	+= $(patsubst %, ${ROOT_ILMBASE}/lib/lib%*.so*, $(LIBS_ILMBASE))

# tbb
ROOT_TBB 		:= ${PIPE_LIBS_ROOT}/tbb/${TBB_VERSION}
LIBS_TBB 		+= tbb
CXXFLAGS 		+= -I${ROOT_TBB}/include
LIBS	 		+= $(patsubst %, -l%, $(LIBS_TBB))
INSTALL_LIBS 	+= $(patsubst %, ${ROOT_TBB}/lib/lib%*.so*, $(LIBS_TBB))

# alembic
ROOT_ALEMBIC	:= ${PIPE_LIBS_ROOT}/alembic/${ALEMBIC_VERSION}
LDFLAGS  		+= -L${ROOT_ALEMBIC}/lib/
CXXFLAGS 		+= -I${ROOT_ALEMBIC}/include
LIBS_ALEMBIC	+= Alembic
INSTALL_LIBS 	+= $(patsubst %, ${ROOT_ALEMBIC}/lib/lib%*.so*, $(LIBS_ALEMBIC))
LIBS_ALEMBIC	+= IECoreAlembic
LIBS	 		+= $(patsubst %, -l%, $(LIBS_ALEMBIC))
INSTALL_LIBS 	+= $(patsubst %, ${ROOT_CORTEX}/alembic/${ALEMBIC_VERSION}/lib/boost${BOOST_VERSION}/lib%*.so, IECoreAlembic)


CXXFLAGS		+= -I/usr/include/x86_64-linux-gnu/
LDFLAGS  		+= -Wl,--enable-new-dtags,-rpath-link=./libs,-rpath=./libs  -L./libs

HEADERS := $(wildcard include/*.h)
SRCS := $(wildcard src/*.cpp)
BIN  := export-isolated-uv-textured split-mesh-using-plane
OBJS := $(patsubst %.cpp, %.o, $(SRCS))
# BINS := $(SRCS:%.cpp=%)

all: docker
all_docker := ${BIN}

libs: libs/.done
libs/.done:
	mkdir -p libs
	cp ${PIPE_LIBS_ROOT}/glew/2.1.0/lib/*.so* ./libs/
	for each in ${INSTALL_LIBS} ; do \
		cp -Lrfuv $$each libs/$$(basename $$each); \
		for dep in $$(ldd $$each | grep atomo | awk '{print $$(NF-1)}') ; do \
			cp -Lrfn $$dep libs/$$(basename $$dep); \
		done ;\
	done ;\
	for each in $$(ls -d ./libs/*) ; do \
		echo $$each ;\
		patchelf --set-rpath ./libs $$each ; \
	done && touch libs/.done

# libs/.post_build_done: ${BIN}
# 	mkdir -p libs/extra
# 	for n in {1..3} ; do \
# 		for each in ${BIN} $$(ls -d ./libs/*) ; do \
# 			for dep in $$(ldd $$each  | grep -v '\./libs' | awk '{print $$(NF-1)}' | egrep -v 'not|ld-linux-x86-64|vdso|libc') ; do \
# 				cp -Lrfu $$dep libs/extra/$$(basename $$dep); \
# 			done ;\
# 		done ;\
# 	done


split-mesh-using-plane: src/split-mesh-using-plane.o
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$$LD_LIBRARY_PATH \
	LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$$LIBRARY_PATH \
	${CXX} -g ${LDFLAGS} ${LIBS} $< -o $@

export-isolated-uv-textured: src/export-isolated-uv-textured.o
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$$LD_LIBRARY_PATH \
	LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$$LIBRARY_PATH \
	${CXX} -g ${LDFLAGS} ${LIBS} $< -o $@

%.o: %.cpp ${HEADERS}
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$$LD_LIBRARY_PATH \
	${CXX} -c ${CXXFLAGS} $< -o $@

_test: ${all_docker}
	./${BIN} tests/2_spheres_test_splitting_texture_map_by_uv*.obj

.PHONY: docker
docker: .PHONY libs
	docker build --compress --force-rm ${DOCKER_ARGS} -f ./docker/Dockerfile -t ubuntu_dev ./docker
	docker run -ti --rm \
		-u $$(id -u ${USER}):$$(id -g ${USER}) \
		-v /etc/passwd:/etc/passwd \
		-v /etc/group:/etc/group \
		-v $$(pwd):/work_dev/ \
		-v /atomo/pipeline/libs/:/atomo/pipeline/libs/ \
	ubuntu_dev /bin/sh -c "cd /work_dev/; make ${all_docker}"

test: docker
	docker run -ti --rm \
		-u $$(id -u ${USER}):$$(id -g ${USER}) \
		-v /etc/passwd:/etc/passwd \
		-v /etc/group:/etc/group \
		-v $$(pwd):/work_dev/ \
		-v /atomo/pipeline/libs/:/atomo/pipeline/libs/ \
	ubuntu_dev /bin/sh -c "cd /work_dev/; make _test"

help:
	@echo ""
	@echo -e "make\t\t- build using docker in ubuntu 20.04"
	@echo -e "make test\t- run tests using the .obj files in tests folder"
	@echo -e "make clean\t- remove everything (compiled code and test tifs), but ./libs folder"
	@echo -e "make nuke\t- remove everything, including ./libs folder"
	@echo -e "make help\t- this help"
	@echo ""

clean:
	@echo "Cleaning up..."
	rm -rvf src/*.o ${BIN} tests/*.tif

nuke: clean
	@echo "Cleaning libraries..."
	rm -rvf  libs
