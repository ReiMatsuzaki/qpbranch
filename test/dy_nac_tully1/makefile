include ../../local.mk
include ${QPBRANCH_ROOT}/conf.mk

-include ${main.d}

main.x: main.o ${libqpbranch}

.PHONY: all clean
all: main.x
	./$<
clean:
	${RM} *.x *.o *.d *.pdf
	${RM} -r _con*
	${RM} -r _fig*

