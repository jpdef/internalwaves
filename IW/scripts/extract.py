# uncompyle6 version 3.7.4
# Python bytecode 3.6 (3379)
# Decompiled from: Python 3.7.7 (default, Mar 26 2020, 15:48:22) 
# [GCC 7.3.0]
# Embedded file name: /home/jpdef/Scripps/Projects/IW/notebooks/scripts/extract.py
# Compiled at: 2020-03-23 14:56:02
# Size of source mod 2**32: 3004 bytes
import tables as pt, pandas as pd, numpy as np, itertools, feather, os, datetime as dt
from tqdm import tqdm
filename = 'llc4320_32-122_20x20x125_soundspeed.mat'

def read_table(file_name):
    return pt.open_file(file_name)


def extract_node(table, node_name):
    full_node_name = '/llc4320/' + node_name
    return table.get_node(full_node_name).read()


def print_table(tab, group=[], nodes=[]):
    for group in tab.walk_groups('/'):
        for array in tab.list_nodes(group, classname='Array'):
            print(array)


def convert_to_dataframe(infile, outpath, tcol, incols, outcols):
    """
    Read in a hierarchial data table and convert to frames
    """
    table = read_table(infile)
    T = extract_node(table, tcol).flatten()
    LAT = extract_node(table, incols[0])[:, 0]
    LON = extract_node(table, incols[1])[0, :]
    DEPTH = extract_node(table, incols[2]).flatten()
    I = list((itertools.product)(*[LAT, LON, DEPTH]))
    zero_padding = int(np.floor(np.log10(len(T))) + 1)
    TC = make_chunks(T, 10)
    for k, tc in progressbar(TC, 'Total'):
        frames = []
        for i, t in enumerate(tc):
            d = {'t': t}
            for j, ic in enumerate(incols):
                d[ic] = [elem[j] for elem in I]

            for oc in outcols:
                O = extract_node(table, oc)
                d[oc] = O[i, :, :].flatten()

            frames.append(pd.DataFrame(d))

        make_featherfiles(frames, 'run', outpath, zp=zero_padding, offset=(k * len(tc)))


def make_featherfiles(frames, bname, dpath, zp=0, offset=0):
    print(offset)
    for t, f in progressbar(frames, 'Writing to Disk'):
        fmt = '{:0>' + str(zp) + '}'
        fname = '%s-%s.fthr' % (bname, fmt.format(t + offset))
        path = os.path.join(dpath, fname)
        feather.write_dataframe(f, path)


def make_chunks(axis, chunk_size):
    timechunks = []
    back_itr = 0
    forward_itr = chunk_size
    while forward_itr < len(axis):
        timechunks.append(axis[back_itr:forward_itr])
        back_itr = forward_itr
        forward_itr += chunk_size

    timechunks.append(axis[back_itr:])
    return timechunks


def progressbar(dataset, desc):
    """
    Desc:
    Helper function that wraps the tqdm library to make 
    function call shorter
    """
    iterator = enumerate(dataset)
    return tqdm(iterator, ascii=True, total=(len(dataset)), leave=True, desc=desc)


def matlab2datetime(matlab_datenum):
    day = dt.datetime.fromordinal(int(matlab_datenum))
    dayfrac = dt.timedelta(days=(matlab_datenum % 1)) - dt.timedelta(days=366)
    return day + dayfrac


if __name__ == '__main__':
    convert_to_dataframe('../matlab/' + filename, '../data/llc', 'time', ['lat', 'lon', 'z'], ['W', 'T', 'S'])
# okay decompiling /home/jpdef/Scripps/Projects/IW/scripts/__pycache__/extract.cpython-36.pyc
