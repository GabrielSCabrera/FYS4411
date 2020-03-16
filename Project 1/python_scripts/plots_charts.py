from file_io import read_energies, read_vals
from matplotlib.ticker import MaxNLocator
from analysis import dataAnalysisClass
import matplotlib.pyplot as plt
import numpy as np
import argparse
import re
import os

def parse_args():

    argparse_desc = ('\033[1mCreates plots and charts based on a dataset.\033[m'
                     '\nExpects a command-line argument \033[3mset\033[m which'
                     ' is a string that points to a dataset directory.  This '
                     'directory should contain 2\033[3mM\033[m files, where '
                     '\033[3mM\033[m is the number of simulations.  The '
                     'contents should consist of \033[3mM\033[m files of the '
                     'format \'\033[3mE_<N>.dat\033[m\', and \033[3mM\033[m '
                     'corresponding \'\033[3mval_<N>.dat\033[m\' files (with '
                     '\033[3m<N>\033[m being the number of \033[3mbosons\033[m '
                     'in the simulation.)')

    set_help = ('path to \033[3mset\033[m of output files.')
    cap_help = ('caption for the generated LaTeX table.')
    fig_help = ('figure number for the generated LaTeX table.')

    parser = argparse.ArgumentParser(description = argparse_desc)

    parser.add_argument('set', type = str, help = set_help)
    parser.add_argument('--cap', type = str, default = None, help = cap_help)
    parser.add_argument('--fig', type = int, default = None, help = fig_help)

    return parser.parse_args()

def prep_plot_energies(path_prefix, N_vals):
    '''
        'path_prefix' should be something like:
        '../data_out/<your_dataset_name>/<your_filename>{:d}.dat'

        Where 'N' in 'N_vals' will replace '{:d}'
    '''
    x_max = 0
    for N in N_vals:
        E = read_energies(path_prefix.format(int(N)))
        if len(E) > x_max:
            x_max = len(E)

    ax = plt.figure().gca()

    x = np.arange(1, x_max+1)
    for m,N in enumerate(N_vals):
        E = read_energies(path_prefix.format(int(N)))
        y = np.ones_like(x)*E[-1]
        cutoff = len(E)
        if int(N) == 1:
            label = f'${int(N):<3d}$ Particle'
        else:
            label = f'${int(N):<3d}$ Particles'
        p = plt.plot(x[:cutoff], E, label = label)
        color = p[0].get_color()
        plt.plot(x[cutoff-1:], y[cutoff-1:], color = color, linestyle = 'dotted')
        if m == len(N_vals)-1:
            plt.plot(x[cutoff-1], y[cutoff-1], 'kx', label = 'Last MC Cycle')
        else:
            plt.plot(x[cutoff-1], y[cutoff-1], 'kx')

    plt.xlabel('Monte-Carlo Cycle')
    plt.ylabel('System Energy')
    plt.grid()
    plt.xlim(1, len(x))
    plt.legend()

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

def fmt_exp_float(x):

    if np.isnan(x):
        return 'NAN'

    string = f'{x:.2E}'

    pattern = r'E([\-\+]?)\d+'
    sign = re.findall(pattern = pattern, string = string)[0]

    pattern = r'E[\-\+]?(\d+)'
    exp = re.findall(pattern = pattern, string = string)[0]
    exp = int(exp)

    pattern = r'(\d+\.\d+)E[\-\+]?\d+'
    val = re.findall(pattern = pattern, string = string)[0]

    if sign == '+':
        sign = ''

    out = f'{val}'
    if exp != 0:
        out += ' \\times 10^{'
        out += sign
        out += f'{exp:d}'
        out += '}'

    return out

def chart_from(path_prefix_E, path_prefix_val, N_vals, caption = None,
               fig_num = None):

    chart = ('\\begin{table}[H]\n'
             '\t\\centering\n'
             '\t\\begin{tabular}{r | c c c}\n'
             '\t\t $N$ Bosons'
             '& $\\left\\langle E \\right\\rangle$ [eV]'
             '& $\\left\\langle \\text{var}(E) \\right\\rangle$ '
             '& $\\left\\langle \\text{std}(E) \\right\\rangle$ \\\\\n'
             '\t\t\\hline\\\\\n')

    for N in N_vals:
        path_E = path_prefix_E.format(int(N))

        print(f'\n\033[1mFILENAME:\033[m {path_E}')

        vals = read_vals(path_prefix_val.format(int(N)))
        anal_dict = get_analyser_data(path_E)

        avgs = []
        vars = []
        stds = []

        for val in anal_dict.values():
            avgs.append(val['avg'])
            vars.append(val['var'])
            stds.append(val['std'])

        avgs = np.array(avgs)
        vars = np.array(vars)
        stds = np.array(stds)

        chart += (f'\t\t${int(N):d}$ & '
                  f'${fmt_exp_float(np.mean(avgs))}$ & '
                  f'${fmt_exp_float(np.mean(vars))}$ & '
                  f'${fmt_exp_float(np.mean(stds))}$\\\\\n')

    chart += '\t\\end{tabular}\n'
    if caption is not None and fig_num is not None:
        chart += '\t\\caption{'
        if caption is not None:
            chart += caption
        if fig_num is not None:
            chart += f'\\label{{table_{fig_num:d}}}'
        chart += '}\n'
    chart += '\\end{table}'

    return chart

def get_analyser_data(path):
    analyzer = dataAnalysisClass(path)
    analyzer.bootstrap()
    analyzer.jackknife()
    analyzer.blocking()
    anal_dict = analyzer.returnOutput()
    return anal_dict

def run_all(cmdline_args):
    all_files = os.listdir(cmdline_args.set)

    if not os.path.isdir(cmdline_args.set + 'plots/'):
        os.mkdir(cmdline_args.set + 'plots/')

    if not os.path.isdir(cmdline_args.set + 'charts/'):
        os.mkdir(cmdline_args.set + 'charts/')

    E_files = []
    val_files = []

    E_N = []
    val_N = []

    E_pat = r'E_\d+\.dat'
    val_pat = r'val_\d+\.dat'

    E_N_pat = r'E_(\d+)\.dat'
    val_N_pat = r'val_(\d+)\.dat'

    for f in all_files:
        if '.dat' in f:
            if re.match(E_pat, f):
                E_files.append(f)
                E_N.append(re.findall(E_N_pat, f)[0])
            elif re.match(val_pat, f):
                val_files.append(f)
                val_N.append(re.findall(val_N_pat, f)[0])
            else:
                raise IOError(f'Invalid datafile \033[3m{f}\033[m found')

    E_N = sorted(list(map(int, E_N)))
    val_N = sorted(list(map(int, val_N)))

    if E_N != val_N or len(E_files) != len(val_files):
        raise IOError('Inconsistent filenames in selected data directory')

    N_vals = E_N
    path_prefix_E = data_dir + 'E_{:d}.dat'

    prep_plot_energies(path_prefix_E, N_vals)
    fig_path = cmdline_args.set + 'plots/'
    plt.savefig(fig_path + 'energies.pdf')

    path_prefix_val = data_dir + 'val_{:d}.dat'
    chart_path = cmdline_args.set + 'charts/'
    chart = chart_from(path_prefix_E, path_prefix_val, N_vals,
                       cmdline_args.set, cmdline_args.fig)
    with open(chart_path + 'tex_table.txt', 'w+') as outfile:
        outfile.write(chart)

if __name__ == '__main__':

    cmdline_args = parse_args()
    data_dir = cmdline_args.set
    msg = f'Invalid Directory \033[3m{data_dir}\033[m'
    assert os.path.isdir(data_dir), msg

    run_all(cmdline_args)
