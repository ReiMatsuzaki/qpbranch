include ../local.mk
include ../conf.mk

GTEST=${EXTERNAL}/googletest/googletest
LIBS += -lgtest -lgtest_main -lpthread
LDFLAGS += -L${GTEST}
INCLUDE += -I${GTEST}/include

CXXFLAGS=-O0 -d -Wall -Wextra

${BUILD}/utest_math.x:  ${BUILD}/mathplus.o  ${BUILD}/utest_math.o
${BUILD}/utest_pwgto.x: ${BUILD}/mathplus.o ${BUILD}/eigenplus.o ${BUILD}/operator.o ${BUILD}/pwgto.o ${BUILD}/utest_pwgto.o
${BUILD}/utest_dy_branch.x: ${BUILD}/mathplus.o ${BUILD}/eigenplus.o ${BUILD}/operator.o ${BUILD}/pwgto.o ${BUILD}/dy_branch.o ${BUILD}/utest_dy_branch.o

.PHONY: check clean
check_%: ${BUILD}/utest_%.x
	$<
#check_math: ${BUILD}/utest_math.x 
#	$<
#check_pwgto:${BUILD}/utest_pwgto.x
#	$<

clean: clean_all

