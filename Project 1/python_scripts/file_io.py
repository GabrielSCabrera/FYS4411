from analysis import dataAnalysisClass
from multiprocessing import Pool
import numpy as np
import re
import os

def get_analyser_data(path):
    analyzer = dataAnalysisClass(path)
    analyzer.runAllAnalyses()
    anal_dict = analyzer.returnOutput()
    return anal_dict

def read_energies(filename):
    pool = Pool()
    with open(filename, 'r') as infile:
        E = infile.read().split()
    E = np.array(pool.map(float, E))
    return E

def read_vals(filename):
    values = {}
    with open(filename, 'r') as infile:
        for line in infile.readlines():
            line = line.split()
            values[line[0]] = float(line[1])
    return values

def read_part_b_file(filename):
    with open(filename, 'r') as infile:
        data = infile.read()
    data = data.split('\n\n\n')
    header = data[0].split('\n')
    runs = data[1].split('\n\n')
    out = []
    headers = {line.split()[0]:int(line.split()[1]) for line in header}
    for n, data in enumerate(runs):
        temp = data.split('\n')
        temp = {line.split()[0]:float(line.split()[1]) for line in temp}
        for key, value in headers.items():
            temp[key] = value
        out.append(temp)
    return out

def read_part_c_file(filename):
    with open(filename, 'r') as infile:
        data = infile.read()
    data = data.split('-------------------------------------')
    header = data[0].split('\n')
    header = list(filter(lambda line: line != "", header))
    runs = []
    for d in range(1, len(data)):
        temp = data[d].split('\n\n')
        for line in temp:
            if line.strip() != '':
                runs.append(line.strip())
    out = []
    headers = {line.split()[0]:int(line.split()[1]) for line in header}
    for n, data in enumerate(runs):
        temp = data.split('\n')
        temp = {line.split()[0]:float(line.split()[1]) for line in temp}
        for key, value in headers.items():
            temp[key] = value
        out.append(temp)
    return out

def read_part_e_file(filename):
    with open(filename, 'r') as infile:
        data = infile.read()
    data = data.split('\n\n\n')
    header = data[0].split('\n')
    runs = data[1].split('\n\n')
    out = []
    headers = {line.split()[0]:int(line.split()[1]) for line in header}
    for n, data in enumerate(runs):
        temp = data.split('\n')
        temp = {line.split()[0]:float(line.split()[1]) for line in temp}
        for key, value in headers.items():
            temp[key] = value
        out.append(temp)
    return out

def read_part_g_file(filename):
    pool = Pool()
    with open(filename, 'r') as infile:
        P = infile.read().split('\n')
    r_range = tuple(re.findall(r'r in \[(\d+), (\d+)\]', P[0])[0])
    P = np.array(pool.map(float, filter(lambda i: i != "", P[1:])))
    out = {'probabilities':P, 'range':(float(r_range[0]), float(r_range[1]))}
    return out

def load_part_b(path):
    all_files = os.listdir(path)
    files = []
    for n,f in enumerate(all_files):
        if '.dat' in f:
            N, dim = re.findall(r'N_(\d+)_dim_(\d+)', f)[0]
            data = read_part_b_file(path + '/' + f)
            for run in data:
                run['N'] = int(N)
                run['dim'] = int(dim)
            files += data
    return files

def load_part_c(path):
    all_files = os.listdir(path)
    files = []
    for n,f in enumerate(all_files):
        if '.dat' in f:
            N, dim = re.findall(r'N_(\d+)_dim_(\d+)', f)[0]
            data = read_part_c_file(path + '/' + f)
            for run in data:
                run['N'] = int(N)
                run['dim'] = int(dim)
            files += data
    return files

def load_part_d(path):
    all_files = os.listdir(path)
    files = []
    prefix = 'alpha'
    N_files = 0
    for f in all_files:
        if prefix in f:
            N_files += 1
    count = 0
    for n,f in enumerate(all_files):
        if prefix in f:
            count += 1
            local_files = os.listdir(path + '/' + f)
            suffixes = []
            N_vals = []
            dim_vals = []
            for l in local_files:
                if '.dat' in l:
                    suffix = re.findall(r'.*(_N_\d+_dt_\d+.dat)', l)[0]
                    N = re.findall(r'.*_N_(\d+)_dt_\d+.dat', l)
                    N = N[0]
                    if suffix not in suffixes:
                        N_vals.append(N)
                        suffixes.append(suffix)
            file_pairs = [[f'E{s}', f'val{s}'] for s in suffixes]
            for m,(i,j) in enumerate(file_pairs):
                out = f'Dir {count:{len(str(N_files))}d}/{N_files} '
                out += f'{m+1:{len(str(len(file_pairs)))}d}/{len(file_pairs)}\t'
                print(out, end = '')
                E = read_energies(path + '/' + f + '/' + i)
                means = get_analyser_data(E)
                val = read_vals(path + '/' + f + '/' + j)
                val['energies'] = E
                val['E'] = means['sample']['avg']
                val['var'] = means['sample']['var']
                val['blocking'] = means['blocking']['avg']
                val['var_blocking'] = means['blocking']['var']
                val['bootstrap'] = means['bootstrap']['avg']
                val['var_bootstrap'] = means['bootstrap']['var']
                val['N'] = int(N_vals[m])
                val['dim'] = 3
                files.append(val.copy())

    return files

def load_part_e(path):
    all_files = os.listdir(path)
    files = []
    for n,f in enumerate(all_files):
        if '.dat' in f:
            N = re.findall(r'N_(\d+)', f)[0]
            data = read_part_e_file(path + '/' + f)
            for run in data:
                run['N'] = int(N)
                run['dim'] = 3
            files += data
    return files

def load_part_g(path):
    all_files = os.listdir(path)
    files = []
    for n,f in enumerate(all_files):
        if '.dat' in f:
            N, potential = re.findall(r'N_(\d+)_(\w+).dat', f)[0]
            data = read_part_g_file(path + '/' + f)
            data['N'] = int(N)
            if potential == 'T':
                data['interacting'] = True
            elif potential == 'OB':
                data['interacting'] = False
            else:
                raise Exception('Unexpected filename format was detected')
            data['dim'] = 3
            files.append(data)
    return files
