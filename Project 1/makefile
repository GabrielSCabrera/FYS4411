CXX = g++
CXXFLAGS = -Wall -g -std=c++17 -fopenmp -O3
# CXXFLAGS = -Wall -g -std=c++17 -openmp -O3

main: main.o ./backend/monte_carlo_class.o ./backend/metropolis.o ./backend/gradient_descent.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o ./matpak/Mat.o #./tests/tests_backend.o
	$(CXX) $(CXXFLAGS) -o main main.o ./matpak/Mat.o ./backend/importance_sampling.o ./backend/monte_carlo_class.o ./backend/metropolis.o ./backend/gradient_descent.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o #./tests/tests_backend.o
	# $(CXX) $(CXXFLAGS) -o main main.cpp ./matpak/Mat.o ./backend/importance_sampling.o ./backend/monte_carlo_class.o ./backend/metropolis.o ./backend/gradient_descent.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o #./tests/tests_backend.o

main.o: main.cpp ./backend/importance_sampling.h ./backend/importance_sampling.o ./backend/monte_carlo_class.h  ./backend/metropolis.h ./backend/gradient_descent.h ./wavefunctions/Psi.h ./wavefunctions/Psi_T.h ./wavefunctions/Psi_OB.h ./matpak/Mat.h #./tests/tests_backend.h
	$(CXX) $(CXXFLAGS) -c main.cpp

./frontend/adagrad.o: ./matpak/Mat.h ./wavefunctions/Psi.h ./wavefunctions/Psi_T.h ./wavefunctions/Psi_OB.h ./backend/monte_carlo.h
./backend/monte_carlo_class.o: ./wavefunctions/Psi.h ./matpak/Mat.h
#./backend/metropolis.o: ./backend/monte_carlo_class.h ./wavefunctions/Psi.h ./matpak/Mat.h
./backend/importance_sampling.o: ./backend/monte_carlo_class.h ./wavefunctions/Psi.h ./matpak/Mat.h
./backend/gradient_decent.o: ./backend/monte_carlo_class.h ./wavefunctions/Psi.h ./matpak/Mat.h
./backend/monte_carlo.o: ./wavefunctions/Psi.h ./matpak/Mat.h
./wavefunctions/Psi.o: ./matpak/Mat.h
./wavefunctions/Psi_OB.o: ./wavefunctions/Psi.h ./matpak/Mat.h
./wavefunctions/Psi_T.o: ./wavefunctions/Psi.h ./matpak/Mat.h
./backend/random.o: ./matpak/Mat.h ./wavefunctions/Psi.h
./matpak/Mat.o: ./matpak/Mat.h
./tests/tests_backend.o: ./tests/tests_backend.h ./wavefunctions/Psi.h ./matpak/Mat.h

valgrind:
	valgrind --leak-check=full ./main

clean:
	rm -f *.o ./frontend/*.o ./backend/*.o ./wavefunctions/*.o ./matpak/*.o #./tests/*.o

test: main
	./main

push: clean main clean
	@ echo "pushing to github"
	@ git add .
	@ git commit -m "Autobackup"
	@ git push origin master

pull: clean
	@ echo "pulling from github"
	@ git pull
