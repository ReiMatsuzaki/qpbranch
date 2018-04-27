#
# see http://www.wagavulin.jp/entry/20130601/1370114999

# -- option --
buildtype:=release

# -- Directories --
SRC=${QPBRANCH_ROOT}/src
BUILD=${QPBRANCH_ROOT}/build/${buildtype}
EXTERNAL=${QPBRANCH_ROOT}/external
EIGEN=${EXTERNAL}/eigen-git-mirror
GTEST=${EXTERNAL}/googletest/googletest
#VPATH=${BUILD}:${SRC}

# -- options --
CXX=c++
INCLUDE=-I${SRC} -I${EIGEN}

ifeq (${buildtype},release)
	CXXFLAGS=-O2 -Wall -Wextra
else ifeq (${buildtype},debug)
	CXXFLAGS += -O0 -g -Wall -Wextran
else
	$(error invalid buildtype)
endif

# -- google test --
GTEST_LIBS = -lpthread
GTEST_CPPFLAGS = -isystem ${GTEST}/include
GTEST_INCLUDE = -I${GTEST}/include -I${GTEST}

# -- compile --
${LIB}:
	mkdir -p $@
${OBJ}:
	mkdir -p $@
${BIN}:
	mkdir -p $@

%.x:
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${CXX} $^ -o $@ ${LIBS} ${LDFLAGS}
%.o: %.cpp ${BUILD}
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $< -o $@
${BUILD}/%.o: ${SRC}/%.cpp
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${CXX} ${CXXFLAGS} ${INCLUDE} -c $< -o $@

${BUILD}/gtest-all.o: ${GTEST}/src/gtest-all.cc
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${CXX} ${CPPFLAGS} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<
${BUILD}/gtest_main.o: ${GTEST}/src/gtest_main.cc
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${CXX} ${CPPFLAGS} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<
${BUILD}/libgtest_main.a: ${BUILD}/gtest-all.o ${BUILD}/gtest_main.o
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${AR} ${ARFLAGS} $@ $^

# -- command --
.PHONY: clean_all 
clean_all:
	${RM} -r ${BUILD}



