# -- Directories --
SRC=${QPBRANCH_ROOT}/src
BUILD=${QPBRANCH_ROOT}/build
EXTERNAL=${QPBRANCH_ROOT}/external
EIGEN=${EXTERNAL}/eigen-git-mirror
VPATH=${BUILD}:${SRC}

# -- options --
CXX=c++
CXXFLAGS=-O2 -Wall
INCLUDE=-I${BUILD} -I${SRC} -I${EIGEN}

# -- utils --
mod2obj = $(addprefix ${BUILD}/, $(addsuffix, .o, $(1)))

# -- compile --

%.x:
	${CXX} $^ -o $@ ${LIBS} ${LDFLAGS}
%.o: %.cpp
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $< -o $@
${BUILD}/%.o: ${SRC}/%.cpp
	@if [ ! -d ${BUILD} ]; \
	   then mkdir -p ${BUILD}; \
	fi
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $< -o $@

# -- command --
.PHONY: clean_all
clean_all:
	${RM} -r ${BUILD}

