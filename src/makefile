include ../local.mk
include ../conf.mk

UTEST=pwgto

${BUILD}/utest_math.x:  $(addprefix mathplus.o utest_math.o)
${BUILD}/utest_pwgto.x: $(addprefix ${BUILD}/, libgtest_main.a mathplus.o eigenplus.o operator.o pwgto.o utest_pwgto.o)
${BUILD}/utest_dy_branch.x: $(addprefix ${BUILD}/, mathplus.o eigenplus.o operator.o pwgto.o dy_branch.o utest_dy_branch.o)



.PHONY: all check

check: ${BUILD}/utest_${UTEST}.x
	$<

all: check

clean: clean_all

