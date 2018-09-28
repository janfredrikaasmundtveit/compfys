# Comment lines
# General makefile for c++ - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CPPflags= c++ -O3
# Here we define the library functions we nee
LIB = -larmadillo -llapack -lblas
# Here we define the name of the executable
PROG= project2.exe
${PROG} :	   	main.o  test-functions.o jacobi.o linalg.o
			${CPPflags} main.o test-functions.o jacobi.o linalg.o ${LIB} -o ${PROG}

main.o :			main.cpp 
		        	${CPPflags} -c main.cpp

test-functitons.o :			test-functions.cpp
		        	${CPPflags} -c test-functions.cpp

jacobi.o :		jacobi.cpp 
		        	${CPPflags} -c jacobi.cpp

linalg.o :		linalg.cpp 
		        	${CPPflags} -c linalg.cpp