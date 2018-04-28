include ../../local.mk
include ${QPBRANCH_ROOT}/conf.mk

-include ${main.d}

OBJS=$(addprefix ${BUILD}/, mathplus.o eigenplus.o operator.o pwgto.o dy_branch.o) main.o

main.o: CPPFLAGS+=${GTEST_CPPFLAGS}
main.o: INCLUDE+=${GTEST_INCLUDE}
main.o: main.cpp

main.x: LIBS+=${GTEST_LIBS}
main.x: ${OBJS}

.PHONY: all

all: main.x
	./$<

clean:
	${RM} *.x *.o *.d

