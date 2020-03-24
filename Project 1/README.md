# Variational Monte Carlo
A collaboration between Gabriel S. Cabrera and Eirill S. Hauge.

This code is for project 1 in the course FYS4411.

The vital parts of the implementation of the Variational Monte Carlo (VMC) is divided into three folders: `matpak`, `variational` and `wavefunctions`.

The files in these folders contain the C++ files for producing raw data for the report. The results are stored in the `results`-folder and further processed using the python scripts in the folder `python_scripts`.

## `matpak`
Contains class representing a matrix used for the positions R.

## `wavefunctions`
The abstract class `Psi` sets up the structure of the classes representing the physical properties of the bosonic systems. Subclass `Psi_OB` represents non-interacting bosons. The sub-class `Psi_T` represents interacting bosons. All calculations of the quantum mechanical properties can be found in these classes. The classes are generalized to handle 1, 2 and 3 dimensions. This makes the classes flexible at the cost of overhead.

These classes use pointers of Mat-objects as parameters to prevent the object from being copied in each function call, saving some CPU time. This may make a significant difference as the methods will be called a large number of times per Monte Carlo simulation.

## `variational`
The abstract superclass `Monte_Carlo` contains the implementation of the Monte Carlo integration. To increase speed, several methods for sampling different properties have been implemented to reduce unnecessary work. The class has two subclasses implementing different methods for the random walk and acceptance ratio.

The subclass `Metropolis` implements the brute-force random walk, while the subclass `Hastings` implements the biased random walk using the quantum force.

The separate function `gradient_descent` optimizes the variational parameter alpha.
