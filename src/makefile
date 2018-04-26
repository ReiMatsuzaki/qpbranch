include ../local.mk
include ../conf.mk

GTEST=${EXTERNAL}/googletest/googletest
LIBS += -lgtest -lgtest_main -lpthread
LDFLAGS += -L${GTEST}
INCLUDE += -I${GTEST}/include

CXXFLAGS=-O0 -d -Wall -Wextra

${BUILD}/utest_pwgto.x: ${BUILD}/mathplus.o ${BUILD}/eigenplus.o ${BUILD}/operator.o ${BUILD}/pwgto.o ${BUILD}/utest_pwgto.o
${BUILD}/utest_math.x:  ${BUILD}/mathplus.o  ${BUILD}/utest_math.o

.PHONY: check clean
check_math: ${BUILD}/utest_math.x 
	$<
check_pwgto:${BUILD}/utest_pwgto.x
	$<

clean: clean_all

