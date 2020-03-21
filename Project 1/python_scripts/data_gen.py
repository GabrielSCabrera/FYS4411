"""
    GENERATE FAKE DATASETS TO TEST PLOTTING FUNCTIONS
"""

import matplotlib.pyplot as plt
import numpy as np
import os

def trend(x, a=1, b=1, c=1, d=0):
    return np.sqrt(c/(a*x+b)) + d

def plot_trend():
    x = np.linspace(0, 1, int(1E3))
    plt.plot(x, trend(x, 1, 3, 1, 0))
    plt.show()

def gen_arrs(N_cycles = 1E3):
    files = { # N PARTICLES                KWARG VALUES
                  '1'        :       {'d' : 0,     'b' : 1},
                  '5'        :       {'d' : 0.01,  'b' : 1.2},
                  '10'       :       {'d' : 0.02,  'b' : 1.4},
                  '50'       :       {'d' : 0.03,  'b' : 1.6},
                  '100'      :       {'d' : 0.04,  'b' : 1.8}
            }

    x = np.linspace(0, 1, int(N_cycles))
    outputs = np.zeros((len(files), len(x)))
    N_vals = []
    for m, (N, kwargs) in enumerate(files.items()):
        t = trend(x, **kwargs)
        rand_interval = np.max(t) - np.min(t)
        rand_interval /= 5
        outputs[m] = t + rand_interval * np.random.random(len(t)) * x**2
        outputs[m] *= 0.5555
        N_vals.append(N)

    return N_vals, outputs

def save_E(N_vals, outputs):
    if not os.path.isdir('../results/test_set/'):
        os.mkdir('../results/test_set/')

    for N, arr in zip(N_vals, outputs):
        fname = f'../results/test_set/E_{int(N):d}.dat'
        np.savetxt(fname = fname, X = arr, fmt = '%f', delimiter = '\n')

def save_val(N_vals, outputs):
    if not os.path.isdir('../results/test_set/'):
        os.mkdir('../results/test_set/')

    N_arr = np.array(list(map(int, N_vals)))
    N_min = np.min(N_arr)
    N_max = np.max(N_arr)
    alpha_0 = 0.5
    beta_0 = 1
    for N, arr in zip(N_vals, outputs):
        fname = f'../results/test_set/val_{int(N):d}.dat'

        alpha = alpha_0-0.02*(float(N)-N_min)/(N_max-N_min)
        beta = beta_0-0.01*(float(N)-N_min)/(N_max-N_min)
        E_avg = np.mean(arr)
        accept = (np.random.random()/50) + 0.95

        outstr = (f'alpha {alpha:f}\n'
                  f'beta {beta:f}\n'
                  f'<E> {E_avg:f}\n'
                  f'accept {accept:f}')

        with open(fname, 'w+') as outfile:
            outfile.write(outstr)

if __name__ == '__main__':
    N_vals, outputs = gen_arrs(1E4)
    save_E(N_vals, outputs)
    save_val(N_vals, outputs)
