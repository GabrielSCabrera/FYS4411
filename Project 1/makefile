CXX = mpic++
CXXFLAGS = -Wall -g -std=c++17

main: main.o ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o ./variational/metropolis.o ./variational/gradient_descent.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o
	#$(CXX) $(CXXFLAGS) -o main main.o ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o ./variational/metropolis.o ./variational/gradient_descent.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o
	$(CXX) $(CXXFLAGS) -o main main.cpp monte_carlo_simulation.cpp ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o ./variational/metropolis.o ./variational/gradient_descent.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o

main.o: main.cpp monte_carlo_simulation.h ./matpak/Mat.h ./variational/hastings.h ./variational/monte_carlo.h  ./variational/metropolis.h ./variational/gradient_descent.h ./wavefunctions/Psi.h ./wavefunctions/Psi_T.h ./wavefunctions/Psi_OB.h
	#$(CXX) $(CXXFLAGS) -c main.cpp 


./variational/monte_carlo.o: ./wavefunctions/Psi.h ./matpak/Mat.h
./variational/metropolis.o: ./variational/monte_carlo.h ./wavefunctions/Psi.h ./matpak/Mat.h
./variational/hastings.o: ./variational/monte_carlo.h ./wavefunctions/Psi.h ./matpak/Mat.h
./variational/gradient_decent.o: ./variational/monte_carlo.h ./wavefunctions/Psi.h ./matpak/Mat.h
./wavefunctions/Psi.o: ./matpak/Mat.h
./wavefunctions/Psi_OB.o: ./wavefunctions/Psi.h ./matpak/Mat.h
./wavefunctions/Psi_T.o: ./wavefunctions/Psi.h ./matpak/Mat.h

./matpak/Mat.o: ./matpak/Mat.h

valgrind:
	valgrind --leak-check=full ./main

clean:
	rm -f *.o ./variational/*.o ./wavefunctions/*.o ./matpak/*.o 

test: main
	make && mpiexec -np 1 ./main


part_b: part_b.cpp ./matpak/Mat.o ./variational/metropolis.o ./variational/monte_carlo.o  ./wavefunctions/Psi_OB.o  ./wavefunctions/Psi.o 
	$(CXX) $(CXXFLAGS) part_b.cpp ./matpak/Mat.o ./variational/metropolis.o ./variational/monte_carlo.o ./wavefunctions/Psi_OB.o ./wavefunctions/Psi.o
	#mpiexec --hostfile hostfile -np 3 ./a.out
	mpiexec -np 2 ./a.out
	make clean && rm a.out

part_c: part_c.cpp ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o  ./wavefunctions/Psi_OB.o  ./wavefunctions/Psi.o 
	$(CXX) $(CXXFLAGS) part_c.cpp ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o ./wavefunctions/Psi_OB.o ./wavefunctions/Psi.o
	#mpiexec --hostfile hostfile -np 3 ./a.out
	mpiexec -np 2 ./a.out
	make clean && rm a.out

part_d: part_d.cpp ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o  ./wavefunctions/Psi_OB.o  ./wavefunctions/Psi.o 
	$(CXX) $(CXXFLAGS) part_d.cpp ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o ./wavefunctions/Psi_OB.o ./wavefunctions/Psi.o
	#mpiexec --hostfile hostfile -np 3 ./a.out
	mpiexec -np 2 ./a.out
	make clean && rm a.out

part_e: part_e.cpp ./matpak/Mat.o ./variational/metropolis.o ./variational/monte_carlo.o  ./wavefunctions/Psi_T.o  ./wavefunctions/Psi.o 
	$(CXX) $(CXXFLAGS) part_e.cpp ./matpak/Mat.o ./variational/metropolis.o ./variational/monte_carlo.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi.o
	mpiexec --hostfile hostfile -np 3 ./a.out
	make clean && rm a.out

part_f: part_f.cpp ./matpak/Mat.o ./variational/metropolis.o ./variational/monte_carlo.o  ./wavefunctions/Psi_T.o  ./wavefunctions/Psi.o ./variational/gradient_descent.o 
	$(CXX) $(CXXFLAGS) part_f.cpp ./matpak/Mat.o ./variational/metropolis.o ./variational/monte_carlo.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi.o  ./variational/gradient_descent.o 
	mpiexec --hostfile hostfile -np 3 ./a.out
	make clean && rm a.out

part_g: part_g.cpp ./matpak/Mat.o ./variational/metropolis.o ./variational/monte_carlo.o  ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o  ./wavefunctions/Psi.o ./variational/gradient_descent.o 
	$(CXX) $(CXXFLAGS) part_g.cpp ./matpak/Mat.o ./variational/metropolis.o ./variational/monte_carlo.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o ./wavefunctions/Psi.o  ./variational/gradient_descent.o 
	mpiexec --hostfile hostfile -np 3 ./a.out
	make clean && rm a.out




push: clean main clean
	@ echo "pushing to github"
	@ git add .
	@ git commit -m "Autobackup"
	@ git push origin master

pull: clean
	@ echo "pulling from github"
	@ git pull
