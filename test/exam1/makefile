include ../../local.mk
include ../../conf.mk

INCLUDE+=-I${SRC}


#OBJS=${BUILD}/mathplus.o ${BUILD}/eigenplus.o ${BUILD}/operator.o ${BUILD}/pwgto.o main.o
OBJS=mathplus.o eigenplus.o operator.o pwgto.o main.o
main.x: ${OBJS}

all: main.x
	./main.x

clean:
	${RM} *.o
	${RM} *.x


