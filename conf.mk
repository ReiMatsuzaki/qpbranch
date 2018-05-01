#
# see http://www.wagavulin.jp/entry/20130601/1370114999

# -- option --
# release or debug
buildtype:=debug
# unit test or calculation name
target:=pwgto

# -- Directories --
SRC:=${QPBRANCH_ROOT}/src
BUILD:=${QPBRANCH_ROOT}/build/${buildtype}
OBJ:=${BUILD}/obj
LIB:=${BUILD}/lib
BIN:=${BUILD}/bin
CALC:=${QPBRANCH_ROOT}/calc
INCLUDE:=${QPBRANCH_ROOT}/include
EXTERNAL:=${QPBRANCH_ROOT}/external
TEST:=${QPBRANCH_ROOT}/test
EIGEN:=${EXTERNAL}/eigen-git-mirror
GTEST:=${EXTERNAL}/googletest/googletest
JSON11:=${EXTERNAL}/json11

# -- dependencies --
SRCS=$(wildcard ${SRC}/*.cpp)
DEPS:=${SRCS:%.cpp=${OBJ}/%.d}
-include ${DEPS}

# -- options --
CXX=c++
CXXFLAGS=-MMD -MP -MF $(@:%.o=%.d)
CPPFLAGS+=-I${INCLUDE} -I${EIGEN} -I${JSON11}

ifeq (${buildtype},release)
  CXXDEBUG=-O2 
  CXXFLAGS=${CXXDEBUG} -Wall -Wextra
else ifeq (${buildtype},debug)
  CXXDEBUG=-O0 -g
  CXXFLAGS+=${CXXDEBUG} -Wall -Wextra
endif

# -- google test --
GTEST_LIBS = -lpthread
GTEST_CPPFLAGS = -isystem ${GTEST}/include -I${GTEST}/include -I${GTEST}

# -- compile --
%.x:
	 @[ -d ${BIN} ] || mkdir -p ${BIN}
	${CXX} $^ -o $@ ${LIBS} ${LDFLAGS}
%.o: %.cpp 
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

${OBJ}/%.o: ${SRC}/%.cpp
	@[ -d ${OBJ} ] || mkdir -p ${OBJ}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

${OBJ}/gtest-all.o: ${GTEST}/src/gtest-all.cc
	@[ -d ${OBJ} ] || mkdir -p ${OBJ}
	${CXX} ${GTEST_CPPFLAGS} ${CXXDEBUG} -c -o $@ $<
${OBJ}/gtest_main.o: ${GTEST}/src/gtest_main.cc
	@[ -d ${OBJ} ] || mkdir -p ${OBJ}
	${CXX} ${GTEST_CPPFLAGS} ${CXXDEBUG} -c -o $@ $<
${LIB}/libgtest_main.a: ${OBJ}/gtest-all.o ${OBJ}/gtest_main.o
	@[ -d ${LIB} ] || mkdir -p ${LIB}
	${AR} ${ARFLAGS} $@ $^

${OBJ}/json11.o: ${JSON11}/json11.cpp
	@[ -d ${OBJ} ] || mkdir -p ${OBJ}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c -o $@ $<

SRCS:=$(wildcard ${SRC}/*.cpp) ${JSON11}/json11.cpp
#OBJS=$(addprefix ${OBJ}/, json11.o mathplus.o eigenplus.o operator.o pwgto_buf.o pwgto1c.o pwgto.o)
OBJS:=$(addprefix ${OBJ}/, $(patsubst %.cpp,%.o,$(notdir ${SRCS})))
${LIB}/libqpbranch.a: ${OBJS}
	@[ -d ${LIB} ] || mkdir -p ${LIB}
	${AR} ${ARFLAGS} $@ $^

# -- command --
.PHONY: check calc clean_all clean_part

check:
	cd ${TEST}/${target}; make all
calc:
	cd ${CALC}/${target}; make all

clean_all:
	${RM} -r ${QPBRANCH_ROOT}/build/

clean_part:
	${RM} -r ${BUILD}



