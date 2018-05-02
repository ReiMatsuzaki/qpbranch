include ../../local.mk
include ${QPBRANCH_ROOT}/conf.mk

-include main.d

main.o: CPPFLAGS+=${GTEST_CPPFLAGS}
main.o: INCLUDE+=${GTEST_INCLUDE}
main.o: main.cpp

main.x: LIBS+=${GTEST_LIBS}
main.x: ${LIB}/libgtest_main.a main.o
	${CXX} main.o -o $@ -lgtest_main ${LIBS} -L${LIB}

.PHONY: all

all: main.x
	./$<

clean:
	${RM} *.x *.o *.d


