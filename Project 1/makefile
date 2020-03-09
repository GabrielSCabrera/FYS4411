CXX = g++
CXXFLAGS = -Wall -g -std=c++17 -fopenmp -O3

main: main.o ./frontend/adagrad.o ./backend/monte_carlo.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o ./backend/random.o ./matpak/Mat.o ./tests/tests_backend.o
	$(CXX) $(CXXFLAGS) -o main main.o ./frontend/adagrad.o ./backend/monte_carlo.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o ./wavefunctions/Psi.o ./backend/random.o ./matpak/Mat.o ./tests/tests_backend.o

main.o: main.cpp ./frontend/adagrad.h ./backend/monte_carlo.h ./wavefunctions/Psi.h ./wavefunctions/Psi_T.h ./wavefunctions/Psi_OB.h ./backend/random.h ./matpak/Mat.h ./tests/tests_backend.h
	$(CXX) $(CXXFLAGS) -c main.cpp

./frontend/adagrad.o: ./matpak/Mat.h ./wavefunctions/Psi.h ./wavefunctions/Psi_T.h ./wavefunctions/Psi_OB.h ./backend/monte_carlo.h
./backend/monte_carlo.o: ./backend/Psi.h ./matpak/Mat.h
./wavefunctions/Psi.o: ./wavefunctions/Psi.h ./matpak/Mat.h
./wavefunctions/Psi_OB.o: ./wavefunctions/Psi.h ./wavefunctions/Psi_OB.h ./matpak/Mat.h
./wavefunctions/Psi_T.o: ./wavefunctions/Psi.h ./wavefunctions/Psi_T.h ./matpak/Mat.h
./backend/random.o: ./matpak/Mat.h ./backend/Psi.h
./matpak/Mat.o: ./matpak/Mat.h
./tests/tests_backend.o: ./tests/tests_backend.h ./backend/Psi.h ./matpak/Mat.h

valgrind:
	valgrind --leak-check=full ./main

clean:
	rm -f *.o ./frontend/*.o ./backend/*.o ./matpak/*.o ./tests/*.o

test: main
	./main
