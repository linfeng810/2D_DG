test: type_declaration.o basis_function.o gmsh_input.o assemble_matrix.o
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fdefault-real-8 tests/test_basis_function.f90 \
		type_declaration.o basis_function.o gmsh_input.o  assemble_matrix.o \
		-o test_sf.exe
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fdefault-real-8 tests/test_gmsh_input.f90 \
		type_declaration.o gmsh_input.o -o test_gmsh.exe

dg_2d: type_declaration.o basis_function.o gmsh_input.o assemble_matrix.o
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fdefault-real-8 main.f90 \
		type_declaration.o basis_function.o gmsh_input.o  assemble_matrix.o \
		-o dg_2d.exe -static

dg_2d_debug: type_declaration.o basis_function.o gmsh_input.o assemble_matrix.o
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fdefault-real-8 main.f90 \
		type_declaration.o basis_function.o gmsh_input.o  assemble_matrix.o \
		-o dg_2d.exe

type_declaration.o: type_declaration.f90 
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fdefault-real-8 -c type_declaration.f90 

basis_function.o: basis_function.f90 
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fdefault-real-8 -c basis_function.f90 

gmsh_input.o: gmsh_input.f90 
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fdefault-real-8 -c gmsh_input.f90 

assemble_matrix.o: assemble_matrix.f90 
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -fdefault-real-8 -c assemble_matrix.f90

clean: 
	rm type_declaration.o basis_function.o gmsh_input.o  assemble_matrix.o *.exe