from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np
import argparse
import file_io
import pandas
import re
import os

label_map = {
             'alpha'    :   '$\\alpha$',
             'E'        :   '$\\left\\langle E \\right\\rangle$',
             'E/Nd'     :   '$\\left\\langle E \\right\\rangle / N d_r$',
             'var'      :   '$\\text{var}(E)$',
             'accept'   :   'Accept Ratio',
             # 'cycles'   :   'Cycles',
             # 'workers'  :   'Workers',
             'N'        :   '$N$',
             'dim'      :   '$d_r$',
             'dt'       :   '$\\Delta t$',
             'blocking' :   '$\\left\\langle E \\right\\rangle_\\text{block}$',
             'var_blocking' :   '$\\left\\langle E \\right\\rangle_\\text{block}$',
             'bootstrap' :   '$\\left\\langle E \\right\\rangle_\\text{boot}$',
             'var_bootstrap' :   '$\\left\\langle E \\right\\rangle_\\text{boot}$',
            }

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

def get_uncertainty(val, var):
    val = fmt_exp_float(val)
    var = fmt_exp_float(var)
    out = f'{f"{val} ± {var}"}'
    out = re.sub(r'±', '\\pm', out)
    return out

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

def fmt_exp_float(x, fmt_f = 'g', fmt_e = '.0E'):

    if isinstance(x, str):
        return x

    if np.isnan(x):
        return 'NAN'

    if int(x) == x and len(str(int(x))) <= 4:
        return str(int(x))
    elif 'e' not in f'{x:{fmt_f}}'.lower():
        return f'{x:{fmt_f}}'

    string = f'{x:{fmt_e}}'

    pattern = r'E([\-\+]?)\d+'
    sign = re.findall(pattern = pattern, string = string)[0]

    pattern = r'E[\-\+]?(\d+)'
    exp = re.findall(pattern = pattern, string = string)[0]
    exp = int(exp)

    pattern = r'(\d+)E[\-\+]?\d+'
    # pattern = r'(\d+\.\d+)E[\-\+]?\d+'
    val = re.findall(pattern = pattern, string = string)[0]

    if sign == '+':
        sign = ''

    if exp != 0:
        out = f'{val}'
        if val != 1:
            out += ' \\times 10^{'
        out += sign
        out += f'{exp:d}'
        out += '}'
    else:
        out = f'{val}'

    return out

def get_col_data(data):
    cols = []
    for d in data:
        for key, val in d.items():
            if key not in cols and key in label_map.keys():
                cols.append(key)

    col_data = {col:[] for col in cols}
    for d in data:
        checked = {col:False for col in cols}
        for key, val in d.items():
            if key in label_map.keys():
                col_data[key].append(val)
                checked[key] = True
        for key, val in checked.items():
            if val is False:
                raise Exception('Inconsistent Labeling in Loaded Files')

    for key, val in col_data.items():
        col_data[key] = np.array(val)

    return col_data

def chart_from_path(path_prefix_E, path_prefix_val, N_vals, caption = None,
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

def sort_col_data(col_data, cols):
    '''
        Sort <dict> 'col_data' by a list of its keys, given by 'col'
    '''
    data = pandas.DataFrame(col_data)
    data = data.sort_values(cols)
    out = {}
    for key in data:
        out[key] = np.array(data[key])
    return out

def divide_by(col_data, col):
    '''
        Split <dict> 'col_data' by one of its keys, given by 'col'
        'col_data' should already be sorted by 'col'
    '''
    data = pandas.DataFrame(col_data)
    sub_data = {}
    for n,i in enumerate(data[col]):
        if i not in sub_data.keys():
            sub_data[i] = [n]
        else:
            sub_data[i].append(n)

    for key, val in sub_data.items():
        sub_data[key] = data.iloc[val].drop(columns = col)
        out = {}

        for d in sub_data[key]:
            out[d] = np.array(sub_data[key][d])
        sub_data[key] = out

    return sub_data

def chart_from_data(col_data, caption = None, fig_num = None, div_col = None):
    '''
        Make sure col_data is sorted by 'div_col' if 'div_col' is not None
    '''
    split_val = None

    chart = ('\\begin{table}[H]\n'
             '\t\\centering\n')

    if caption is not None and fig_num is not None:
        chart += '\t\\caption{'
        if caption is not None:
            chart += caption
        if fig_num is not None:
            chart += f'\\label{{table_{fig_num}}}'
        chart += '}\n'

    chart += '\t\\begin{tabular}{'

    cols = list(col_data.keys())

    if div_col is not None:
        chart += 'c | ' + 'c '*(len(cols)-1)
    else:
        chart += 'c '*len(cols)

    chart = chart[:-1]
    chart += '}\n\t\t'
    if div_col is not None:
        chart += f'{label_map[div_col]} & '
    for col in cols:
        if div_col is None or col != div_col:
            chart += f'{label_map[col]} & '
    chart = chart[:-2] + '\\\\\n\t\t\\hline\n\t\t\\hline\n'

    for i in range(len(col_data[cols[0]])):
        chart += '\t\t'
        if div_col is not None:
            if split_val is None:
                split_val = col_data[div_col][0]
            else:
                if split_val != col_data[div_col][i]:
                    split_val = col_data[div_col][i]
                    chart += '\\hline\n\t\t'
            chart += f'${fmt_exp_float(col_data[div_col][i])}$ & '

        for col in cols:
            if col == 'alpha':
                chart += f'${fmt_exp_float(col_data[col][i], ".4f")}$ & '
            elif col == 'accept':
                if f'{col_data[col][i]:.2f}' == f'{1:.2f}':
                    chart += f'$0.99$ & '
                else:
                    chart += f'${fmt_exp_float(col_data[col][i], ".2f")}$ & '
            elif col == 'dt':
                chart += f'${fmt_exp_float(col_data[col][i], ".3f")}$ & '
            elif div_col is None or col != div_col:
                chart += f'${fmt_exp_float(col_data[col][i])}$ & '

        chart = chart[:-2] + '\\\\\n'

    chart += '\t\\end{tabular}\n'
    chart += '\\end{table}'

    return chart

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

def combine_cols_uncertainty(data, vals_label, vars_label):
    combined = []
    dataframe = pandas.DataFrame(data)
    for val, var in zip(dataframe[vals_label], dataframe[vars_label]):
        combined.append(get_uncertainty(val, var))
    col_vals = []
    for i in range(len(data[vals_label])):
        col_vals.append(combined[i])
    col_vals = np.array(col_vals)
    data[vals_label] = col_vals
    return data

def part_b(path, main_dir, show = False):

    data = file_io.load_part_b(path)
    data = get_col_data(data)
    if show:
        for run in data:
            for key, val in run.items():
                print(f'{key:10s}\t{val:g}')
            print()

    data['std_1'] = np.sqrt(data['var'])
    data = combine_cols_uncertainty(data, 'E', 'std_1')
    del data['std_1']
    data['std_2'] = np.sqrt(data['var']/(data['N']*data['dim']))
    data = combine_cols_uncertainty(data, 'E/Nd', 'std_2')
    del data['std_2']
    del data['var']

    data = sort_col_data(data, ['dim', 'N', 'alpha'])
    sub_data = divide_by(data, 'dim')

    out = ''
    for n, (key, data) in enumerate(sub_data.items()):
        div_col = 'N'
        caption = f'PART B, {label_map["dim"]}={key}'
        fig_num = '1.' + str(n)
        chart = chart_from_data(data, caption, fig_num, div_col)
        out += chart + '\n'
    with open(path + '/charts/tex_table.txt', 'w+') as outfile:
        outfile.write(out)
    with open(main_dir + '/all_tables.txt', 'w+') as outfile:
        outfile.write(out)

def part_c(path, main_dir, show = False):
    data = file_io.load_part_c(path)
    data = get_col_data(data)
    if show:
        for run in data:
            for key, val in run.items():
                print(f'{key:10s}\t{val:g}')
            print()

    data['std_1'] = np.sqrt(data['var'])
    data = combine_cols_uncertainty(data, 'E', 'std_1')
    del data['std_1']
    data['std_2'] = np.sqrt(data['var']/(data['N']*data['dim']))
    data = combine_cols_uncertainty(data, 'E/Nd', 'std_2')
    del data['std_2']
    del data['var']

    data = sort_col_data(data, ['dim','N','alpha','dt'])
    sub_data = divide_by(data, 'dim')

    out = ''
    for n, (key, data) in enumerate(sub_data.items()):
        div_col = 'N'
        caption = f'PART C, {label_map["dim"]}={key}'
        fig_num = '2.' + str(n)
        chart = chart_from_data(data, caption, fig_num, div_col)
        out += chart + '\n'
    with open(path + '/charts/tex_table.txt', 'w+') as outfile:
        outfile.write(out)
    with open(main_dir + '/all_tables.txt', 'a+') as outfile:
        outfile.write(out)

def part_d(path, main_dir, show = False):
    data = file_io.load_part_d(path)
    out = ''

    data = get_col_data(data)
    if show:
        for run in data:
            for key, val in run.items():
                if isinstance(val, np.ndarray):
                    print(f'{key:10s}\t{"array"}')
                else:
                    print(f'{key:10s}\t{val:g}')
            print()

    data = sort_col_data(data, ['N','alpha','dt'])

    data['std'] = np.sqrt(data['var_blocking'])
    data = combine_cols_uncertainty(data, 'blocking', 'std')

    data['std'] = np.sqrt(data['var_bootstrap'])
    data = combine_cols_uncertainty(data, 'bootstrap', 'std')
    del data['std']

    del data['E']
    del data['var']
    del data['var_blocking']
    del data['var_bootstrap']
    del data['dim']

    div_col = 'N'
    caption = f'PART D'
    fig_num = 3
    chart = chart_from_data(data, caption, fig_num, div_col)
    out += chart + '\n'
    with open(path + '/charts/tex_table.txt', 'w+') as outfile:
        outfile.write(out)
    with open(main_dir + '/all_tables.txt', 'a+') as outfile:
        outfile.write(out)

def part_e(path, main_dir, show = False):
    data = file_io.load_part_e(path)
    data = get_col_data(data)
    if show:
        for run in data:
            for key, val in run.items():
                if isinstance(val, np.ndarray):
                    print(f'{key:10s}\t{val}')
                else:
                    print(f'{key:10s}\t{val:g}')
            print()

    data['std_1'] = np.sqrt(data['var'])
    data = combine_cols_uncertainty(data, 'E', 'std_1')
    del data['std_1']
    data['std_2'] = np.sqrt(data['var']/(data['N']*data['dim']))
    data = combine_cols_uncertainty(data, 'E/Nd', 'std_2')
    del data['std_2']
    del data['var']

    data = sort_col_data(data, ['dim','N','alpha'])
    sub_data = divide_by(data, 'dim')

    out = ''
    for n, (key, data) in enumerate(sub_data.items()):
        div_col = 'N'
        caption = f'PART E, {label_map["dim"]}={key}'
        fig_num = '4.' + str(n)
        chart = chart_from_data(data, caption, fig_num, div_col)
        out += chart + '\n'
    with open(path + '/charts/tex_table.txt', 'w+') as outfile:
        outfile.write(out)
    with open(main_dir + '/all_tables.txt', 'a+') as outfile:
        outfile.write(out)

def part_f(path, main_dir, show = False):
    data = file_io.load_part_e(path)
    data = get_col_data(data)
    if show:
        for run in data:
            for key, val in run.items():
                if isinstance(val, np.ndarray):
                    print(f'{key:10s}\t{val}')
                else:
                    print(f'{key:10s}\t{val:g}')
            print()

    data['std_1'] = np.sqrt(data['var'])
    data = combine_cols_uncertainty(data, 'E', 'std_1')
    del data['std_1']
    data['std_2'] = np.sqrt(data['var']/(data['N']*data['dim']))
    data = combine_cols_uncertainty(data, 'E/Nd', 'std_2')
    del data['std_2']
    del data['var']
    del data['dim']
    
    data = sort_col_data(data, ['N','alpha'])

    out = ''
    for n, (key, data) in enumerate(sub_data.items()):
        div_col = 'N'
        caption = f'PART F, {label_map["dim"]}={key}'
        fig_num = '5.' + str(n)
        chart = chart_from_data(data, caption, fig_num, div_col)
        out += chart + '\n'
    with open(path + '/charts/tex_table.txt', 'w+') as outfile:
        outfile.write(out)
    with open(main_dir + '/all_tables.txt', 'a+') as outfile:
        outfile.write(out)

def part_g(path, main_dir, show = False):
    data = file_io.load_part_g(path)
    data = get_col_data(data)
    if show:
        for run in data:
            for key, val in run.items():
                if not isinstance(val, (float, int)):
                    print(f'{key:12s}\t{val}')
                else:
                    print(f'{key:12s}\t{val:g}')
            print()
    chart = chart_from_data(data, caption = 'Part G', fig_num = 6)
    with open(path + '/charts/tex_table.txt', 'w+') as outfile:
        outfile.write(chart)

def run_all_parts(cmdline_args, show = False):
    all_files = os.listdir(cmdline_args.set)
    parts = {
             'part_b'   :   part_b,
             'part_c'   :   part_c,
             'part_d'   :   part_d,
             'part_e'   :   part_e,
             'part_f'   :   part_f,
             # 'part_g'   :   part_g
            }

    for part, func in parts.items():
        if part not in all_files:
            print(f'Missing directory \'{part}\', skipping to next part...')
            continue
        if show:
            msg = ' '.join(part.upper().split('_'))
            print(f'\033[1m{msg}\033[m START')

        path = cmdline_args.set + part

        if not os.path.isdir(path + '/charts/'):
            os.mkdir(path + '/charts/')

        if not os.path.isdir(path + '/plots/'):
            os.mkdir(path + '/plots/')

        parts[part](path, cmdline_args.set, show)
        if show:
            print(f'\033[1m{msg}\033[m FINISH')

if __name__ == '__main__':

    cmdline_args = parse_args()
    data_dir = cmdline_args.set
    msg = f'Invalid Directory \033[3m{data_dir}\033[m'
    assert os.path.isdir(data_dir), msg

    run_all_parts(cmdline_args, show = False)
