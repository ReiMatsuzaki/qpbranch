include ../../local.mk
include ${QPBRANCH_ROOT}/conf.mk

-include ${main.d}

main.o: CPPFLAGS+=${GTEST_CPPFLAGS}
main.o: INCLUDE+=${GTEST_INCLUDE}
main.x: LIBS+=${GTEST_LIBS}
main.x: main.o ${libqpbranch} ${libgtest_main}

.PHONY: all clean plot
all: main.x
	./$<

plot:
	python plot.py

clean:
	${RM} *.x *.o *.d
	${RM} -r _con_*
	${RM} -r _fig_*

