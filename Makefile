test: type_declaration.o basis_function.o gmsh_input.o
	gfortran tests/test_basis_function.f90 \
		type_declaration.o basis_function.o gmsh_input.o -o a.out

type_declaration.o: type_declaration.f90 
	gfortran -c type_declaration.f90 

basis_function.o: basis_function.f90 
	gfortran -c basis_function.f90 

gmsh_input.o: gmsh_input.f90 
	gfortran -c gmsh_input.f90 

clean: 
	rm type_declaration.o basis_function.o gmsh_input.o