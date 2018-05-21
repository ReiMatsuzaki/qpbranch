include ../../local.mk
include ${QPBRANCH_ROOT}/conf.mk

-include ${main.d}

main.o: CPPFLAGS+=${GTEST_CPPFLAGS}
main.o: INCLUDE+=${GTEST_INCLUDE}
main.o: main.cpp

main.x: LIBS+=${GTEST_LIBS}
main.x: ${LIB}/libqpbranch.a ${LIB}/libgtest_main.a main.o
	${CXX} main.o -o $@ -L${LIB} -lqpbranch -lgtest_main ${LIBS} ${LDFLAGS}

.PHONY: all

all: main.x
	./$<

clean:
	${RM} *.x *.o *.d


