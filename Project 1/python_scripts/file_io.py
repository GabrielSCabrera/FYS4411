from multiprocessing import Pool
import numpy as np

def read_energies(filename):
    pool = Pool()
    with open(filename, 'r') as infile:
        E = infile.read().split()
    E = np.array(pool.map(float, E))
    return E

def read_vals(filename):
    with open(filename, 'r') as infile:
        values = {}
        for line in infile.readlines():
            line = line.split()
            values[line[0]] = float(line[1])
    return values

if __name__ == '__main__':

    """LOADING TEST FILES"""

    data_dir = '../data_out/'
    test_E_file = 'test_E.dat'
    test_val_file = 'test_val.dat'

    E = read_energies(data_dir + test_E_file)
    values = read_vals(data_dir + test_val_file)

    print(E)
    print(values)
