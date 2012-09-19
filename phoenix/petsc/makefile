# tested with petsc3.2-p5

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

ex9: ex9.o chkopts
	-${CLINKER} -o ex9 ex9.o ${PETSC_SNES_LIB}
	${RM} ex9.o

ex9-petsc3p2: ex9-petsc3p2.o chkopts
	-${CLINKER} -o ex9-petsc3p2 ex9-petsc3p2.o ${PETSC_SNES_LIB}
	${RM} ex9-petsc3p2.o

obstacle: obstacle.o chkopts
	-${CLINKER} -o obstacle obstacle.o ${PETSC_SNES_LIB}
	${RM} obstacle.o

wpi: wpi.o chkopts
	-${CLINKER} -o wpi wpi.o ${PETSC_SNES_LIB}
	${RM} wpi.o


