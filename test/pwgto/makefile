include ../../local.mk
include ${QPBRANCH_ROOT}/conf.mk

-include ${main.d}

main.o: CPPFLAGS+=${GTEST_CPPFLAGS}
main.o: INCLUDE+=${GTEST_INCLUDE}
main.x: LIBS+=${GTEST_LIBS}
main.x: main.o ${libqpbranch} ${libgtest_main}

.PHONY: all clean
all: main.x
	./main.x
clean:
	${RM} *.x *.o *.d


