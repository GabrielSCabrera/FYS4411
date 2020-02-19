CXX = g++
CXXFLAGS = -Wall -g -std=c++17 -fopenmp -O3

main: main.o ./frontend/project_1b.o ./backend/Psi.o ./backend/random.o ./matpak/Mat.o ./matpak/tools.o ./matpak/Vec.o ./tests/tests_backend.o ./tests/tests_matpak.o
	$(CXX) $(CXXFLAGS) -o main main.o ./frontend/project_1b.o ./backend/Psi.o ./backend/random.o ./matpak/Mat.o ./matpak/tools.o ./matpak/Vec.o ./tests/tests_backend.o ./tests/tests_matpak.o

main.o: main.cpp ./frontend/project_1b.h ./backend/Psi.h ./backend/random.h ./matpak/Mat.h ./matpak/tools.h ./matpak/Vec.h ./tests/tests_backend.h ./tests/tests_matpak.h
	$(CXX) $(CXXFLAGS) -c main.cpp

./frontend/project_1b.o: ./backend/Psi.h ./matpak/Mat.h
./backend/Psi.o: ./backend/Psi.h ./matpak/Mat.h
./backend/random.o: ./matpak/Mat.h
./matpak/Mat.o: ./matpak/Mat.h ./matpak/Vec.h ./matpak/tools.h
./matpak/tools.o: ./matpak/tools.h ./matpak/Mat.h ./matpak/Vec.h
./matpak/Vec.o: ./matpak/Vec.h ./matpak/tools.h
./tests/tests_backend.o: ./tests/tests_backend.h ./backend/Psi.h ./matpak/Mat.h
./tests/tests_matpak.o: ./tests/tests_matpak.h ./matpak/Mat.h ./matpak/Vec.h ./matpak/tools.h

valgrind:
	valgrind --leak-check=full ./main

clean:
	rm -f *.o ./frontend/*.o ./backend/*.o ./matpak/*.o ./tests/*.o
