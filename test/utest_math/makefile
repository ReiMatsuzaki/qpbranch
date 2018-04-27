include ../../local.mk
include ${QPBRANCH_ROOT}/conf.mk

-include ${main.d}

OBJS=$(addprefix ${BUILD}/, libgtest_main.a mathplus.o) main.o

main.o: CPPFLAGS+=${GTEST_CPPFLAGS}
main.o: main.cpp

main.x: LIBS+=${GTEST_LIBS}
main.x: ${BUILD}/libgtest_main.a ${BUILD}/mathplus.o main.o

.PHONY: all

all: main.x
	./$<

clean:
	${RM} *.x *.o *.d

