
CFLAGS	        =
FFLAGS	        =
CPPFLAGS        =
FPPFLAGS        =
LOCDIR          = 
EXAMPLESC       = Coagulation.c
EXAMPLESF       = 
EXAMPLESFH      = 
MANSEC          = 

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

coagulation: Coagulation.o chkopts
	-${CLINKER} -g -o Coagulation Coagulation.o ${PETSC_TS_LIB}
	${RM} Coagulation.o

coagulationo3: Coagulation.o chkopts
	-${CLINKER} -o Coagulation Coagulation.o ${PETSC_TS_LIB} -O3
	${RM} -f Coagulation.o

#--------------------------------------------------------------------------------

clean::
	${RM} Coagulation

#---------------------------------------------------------------------------------

runcoagulation:
	-@${MPIEXEC} -n 4 ./Coagulation > Coagulation.tmp 2>&1; \

runcoagulation_withlog:
	-@${MPIEXEC} -n 4 ./Coagulation -log_summary > Coagulation.tmp 2>&1; \
