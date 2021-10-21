test: type_declaration.o basis_function.o gmsh_input.o assemble_matrix.o dirlnsol.o
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow tests/test_basis_function.f90 \
		type_declaration.o basis_function.o gmsh_input.o  assemble_matrix.o \
		-o test_sf.exe
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow tests/test_gmsh_input.f90 \
		type_declaration.o gmsh_input.o -o test_gmsh.exe
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow tests/test_solver.f90 \
		dirlnsol.o -o test_solver.exe

dg_2d: type_declaration.o basis_function.o gmsh_input.o assemble_matrix.o dirlnsol.o
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow main.f90 \
		type_declaration.o basis_function.o gmsh_input.o  assemble_matrix.o dirlnsol.o\
		-o dg_2d.exe -static

dg_2d_debug: type_declaration.o basis_function.o gmsh_input.o assemble_matrix.o dirlnsol.o
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow main.f90 \
		type_declaration.o basis_function.o gmsh_input.o  assemble_matrix.o dirlnsol.o\
		-o dg_2d.exe

type_declaration.o: type_declaration.f90 
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -c type_declaration.f90 

basis_function.o: basis_function.f90 
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -c basis_function.f90 

gmsh_input.o: gmsh_input.f90 
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -c gmsh_input.f90 

assemble_matrix.o: assemble_matrix.f90 
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -c assemble_matrix.f90

dirlnsol.o: dirlnsol.f90
	gfortran -g -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -c dirlnsol.f90

clean: 
	@rm *.o *.exe *.mod
