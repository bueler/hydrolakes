# tested with petsc3.2-p5

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

obstacle: obstacle.o chkopts
	-${CLINKER} -o obstacle obstacle.o ${PETSC_SNES_LIB}
	${RM} obstacle.o

wpi: wpi.o chkopts
	-${CLINKER} -o wpi wpi.o ${PETSC_SNES_LIB}
	${RM} wpi.o


