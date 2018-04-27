#
# see http://www.wagavulin.jp/entry/20130601/1370114999

# -- option --
# release or debug
buildtype:=debug
# unit test name
utest:=pwgto

# -- Directories --
SRC:=${QPBRANCH_ROOT}/src
BUILD:=${QPBRANCH_ROOT}/build/${buildtype}
INCLUDE:=${QPBRANCH_ROOT}/include
EXTERNAL:=${QPBRANCH_ROOT}/external
TEST:=${QPBRANCH_ROOT}/test
EIGEN:=${EXTERNAL}/eigen-git-mirror
GTEST:=${EXTERNAL}/googletest/googletest
JSON11:=${EXTERNAL}/json11
#VPATH=${BUILD}:${SRC}
#VPATH:=${SRC}

# -- dependencies --
SRCS=$(wildcard ${SRC}/*.cpp)
DEPS:=${SRCS:%.cpp=${BUILD}/%.d}
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
	 @[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${CXX} $^ -o $@ ${LIBS} ${LDFLAGS}
%.o: %.cpp 
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

${BUILD}/%.o: ${SRC}/%.cpp
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

${BUILD}/gtest-all.o: ${GTEST}/src/gtest-all.cc
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${CXX} ${GTEST_CPPFLAGS} ${CXXDEBUG} -c -o $@ $<
${BUILD}/gtest_main.o: ${GTEST}/src/gtest_main.cc
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${CXX} ${GTEST_CPPFLAGS} ${CXXDEBUG} -c -o $@ $<
${BUILD}/libgtest_main.a: ${BUILD}/gtest-all.o ${BUILD}/gtest_main.o
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${AR} ${ARFLAGS} $@ $^

${BUILD}/json11.o: ${JSON11}/json11.cpp
	@[ -d ${BUILD} ] || mkdir -p ${BUILD}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c -o $@ $<

# -- command --
.PHONY: clean_all check

check:
	cd ${TEST}/utest_${utest}; make all

clean_all:
	${RM} -r ${QPBRANCH_ROOT}/build

clean_part:
	${RM} -r ${BUILD}



