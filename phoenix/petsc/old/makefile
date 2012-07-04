# tested with petsc3.2-p5

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

onestep: onestep.o chkopts
	-${CLINKER} -o onestep onestep.o ${PETSC_SNES_LIB}
	${RM} onestep.o

porous: porous.o chkopts
	-${CLINKER} -o porous porous.o ${PETSC_TS_LIB}
	${RM} porous.o

bh: bh.o chkopts
	-${CLINKER} -o bh bh.o ${PETSC_TS_LIB}
	${RM} bh.o

kw: kw.o chkopts
	-${CLINKER} -o kw kw.o ${PETSC_TS_LIB}
	${RM} kw.o

wph: wph.o chkopts
	-${CLINKER} -o wph wph.o ${PETSC_TS_LIB}
	${RM} wph.o

