include ../local.mk
include ../conf.mk

GTEST=${EXTERNAL}/googletest/googletest
LIBS += -lgtest -lgtest_main -lpthread
LDFLAGS += -L${GTEST}
INCLUDE += -I${GTEST}/include

${BUILD}/utest_pwgto.x: ${BUILD}/pwgto.o ${BUILD}/utest_pwgto.o
${BUILD}/utest_math.x:  ${BUILD}/mathplus.o  ${BUILD}/utest_math.o

.PHONY: check clean
check: ${BUILD}/utest_math.x ${BUILD}/utest_pwgto.x 
	$<

clean:
	${RM} *.x
	${RM} *.o
