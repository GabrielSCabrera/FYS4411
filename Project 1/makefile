CXX = g++
CXXFLAGS = -Wall -g -std=c++17

main: main.o ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o ./variational/metropolis.o ./variational/gradient_descent.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o
	#$(CXX) $(CXXFLAGS) -o main main.o ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o ./variational/metropolis.o ./variational/gradient_descent.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o
	$(CXX) $(CXXFLAGS) -o main main.cpp ./matpak/Mat.o ./variational/hastings.o ./variational/monte_carlo.o ./variational/metropolis.o ./variational/gradient_descent.o ./wavefunctions/Psi.o ./wavefunctions/Psi_T.o ./wavefunctions/Psi_OB.o

main.o: main.cpp  ./matpak/Mat.h ./variational/hastings.h ./variational/monte_carlo.h  ./variational/metropolis.h ./variational/gradient_descent.h ./wavefunctions/Psi.h ./wavefunctions/Psi_T.h ./wavefunctions/Psi_OB.h
	$(CXX) $(CXXFLAGS) -c main.cpp

#./frontend/adagrad.o: ./matpak/Mat.h ./wavefunctions/Psi.h ./wavefunctions/Psi_T.h ./wavefunctions/Psi_OB.h ./backend/monte_carlo.h
./variational/monte_carlo.o: ./wavefunctions/Psi.h ./matpak/Mat.h
./variational/metropolis.o: ./variational/monte_carlo.h ./wavefunctions/Psi.h ./matpak/Mat.h
./variational/hastings.o: ./variational/monte_carlo.h ./wavefunctions/Psi.h ./matpak/Mat.h
./variational/gradient_decent.o: ./variational/monte_carlo.h ./wavefunctions/Psi.h ./matpak/Mat.h
#./backend/monte_carlo.o: ./wavefunctions/Psi.h ./matpak/Mat.h
./wavefunctions/Psi.o: ./matpak/Mat.h
./wavefunctions/Psi_OB.o: ./wavefunctions/Psi.h ./matpak/Mat.h
./wavefunctions/Psi_T.o: ./wavefunctions/Psi.h ./matpak/Mat.h
#./backend/random.o: ./matpak/Mat.h ./wavefunctions/Psi.h
./matpak/Mat.o: ./matpak/Mat.h
#./tests/tests_backend.o: ./tests/tests_backend.h ./wavefunctions/Psi.h ./matpak/Mat.h

valgrind:
	valgrind --leak-check=full ./main

clean:
	rm -f *.o ./frontend/*.o ./variational/*.o ./wavefunctions/*.o ./matpak/*.o 

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
