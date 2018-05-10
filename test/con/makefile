include ../../local.mk
include ${QPBRANCH_ROOT}/conf.mk

-include ${main.d}

main.o: CPPFLAGS+=${GTEST_CPPFLAGS}
main.o: INCLUDE+=${GTEST_INCLUDE}
main.o: main.cpp
main.x: LIBS+=${GTEST_LIBS}
main.x: main.o

.PHONY: all clean
all: libqpbranch libgtest_main main.x
	${RM} -r _con
	./main.x
	cat _con/z/0001.json

clean:
	${RM} *.x *.o *.d
	${RM} -r _con

