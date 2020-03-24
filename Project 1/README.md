# Variational Monte Carlo 
A collaboration between Gabriel S. Goater and Eirill S. Hauge.


This code is for project 1 in the course FYS4411. 

The vital parts of the implementation of variational Monte Carlo is divided into three folders: `matpak`, `variational` and `wavefunctions`. 

The files in the this folder contain the c++-files for producing raw data for the report. The results are stored in the `result`-folder and further processed using the python scripts in the folder `python_scripts`. 

## `matpack` 
Contains class representing a matrix used for the positions R.

## `wavefunctions` 
The super class `Psi` sets up the structure of the classes representing the physical properties of the bosonic systems. Sub-class `Psi_OB` represents non-interacting bosons. The sub-class `Psi_T` represents interacting bosons. All calculations of the quantum mechanical properties can be found in these classes. The classes are generalized to handle 1, 2 and 3 dimensions. This makes the classes flexible, at the cost of speed.

These classes use pointers of Mat-objects as parameters to prevent the object from being copied in each function call, saving some CPU time. This may make a sigificant difference as the methods will be called a large number of times per Monte Carlo simulation.

## `variational` 
The abstract super class `Monte_Carlo` containes the implementation of the Monte Carlo integration. To increase speed, several methods for sampling different properties have been implemented to reduce unnecessary work. The class has two subclasses implementing different methods for the random walk and acceptance ratio.

The sub-class `Metropolis` implements the brute force random walk, while the sub-class `Hastings` implements the biased random walk using the quantum force. 

The seperate function `gradient_decent` optimizes the variational parameter alpha.