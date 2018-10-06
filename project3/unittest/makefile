# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG= test.exe
${PROG} :	   	test-main.o  test-functions.o  factorial.o 
			${CPPflags} test-main.o test-functions.o factorial.o ${LIB} -o ${PROG}



test-main.o :			test-main.cpp
		        	${CPPflags} -c test-main.cpp

test-functitons.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp


