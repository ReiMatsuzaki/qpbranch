include ../../local.mk
include ${QPBRANCH_ROOT}/conf.mk

-include ${main.d}

#OBJS=$(addprefix ${BUILD}/, libgtest_main.a json11.o mathplus.o eigenplus.o operator.o pwgto.o pwgto_buf.o pwgto1c.o) main.o

main.o: CPPFLAGS+=${GTEST_CPPFLAGS}
main.o: INCLUDE+=${GTEST_INCLUDE}
main.o: main.cpp

main.x: LIBS+=${GTEST_LIBS}
main.x: ${BUILD}/libqpbranch.a ${BUILD}/libgtest_main.a main.o
	${CXX} main.o -o $@ -L${BUILD} -lqpbranch -lgtest_main  ${LIBS} ${LDFLAGS}

.PHONY: all

all: main.x
	./$<

clean:
	${RM} *.x *.o *.d


