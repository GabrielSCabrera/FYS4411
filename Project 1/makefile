CXX = g++
CXXFLAGS = -Wall -g -std=c++17 -fopenmp -O3

main: main.o ./frontend/monte_carlo.o ./backend/Psi.o ./backend/random.o ./matpak/Mat.o ./tests/tests_backend.o
	$(CXX) $(CXXFLAGS) -o main main.o ./frontend/monte_carlo.o ./backend/Psi.o ./backend/random.o ./matpak/Mat.o ./tests/tests_backend.o

main.o: main.cpp ./frontend/monte_carlo.h ./backend/Psi.h ./backend/random.h ./matpak/Mat.h ./tests/tests_backend.h
	$(CXX) $(CXXFLAGS) -c main.cpp

./frontend/monte_carlo.o: ./backend/Psi.h ./matpak/Mat.h
./backend/Psi.o: ./backend/Psi.h ./matpak/Mat.h
./backend/random.o: ./matpak/Mat.h ./backend/Psi.h
./matpak/Mat.o: ./matpak/Mat.h
./tests/tests_backend.o: ./tests/tests_backend.h ./backend/Psi.h ./matpak/Mat.h

valgrind:
	valgrind --leak-check=full ./main

clean:
	rm -f *.o ./frontend/*.o ./backend/*.o ./matpak/*.o ./tests/*.o

test: main
	./main
