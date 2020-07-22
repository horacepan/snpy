import time
import pdb
import copy

class YoungTableau:
    def __init__(self, yt_lst, n):
        self._yt_lst = yt_lst
        self._rows = {}
        self._cols = {}
        self._fill_rows_cols()

    def _fill_rows_cols(self):
        for ridx, row in enumerate(self._yt_lst):
            for cidx, val in enumerate(row):
                self._rows[val] = ridx
                self._cols[val] = cidx

    def get_row(self, x):
        return self._rows[x]

    def get_col(self, x):
        return self._cols[x]

    def ax_dist(self, i, j):
        ri = self.get_row(i)
        rj = self.get_row(j)
        ci = self.get_row(i)
        cj = self.get_row(j)

        return (cj - ci) - (rj - ri)

def gen_standard_tableaux(partition):
    n = sum(partition)
    tabs = [[]]

    for i in range(n):
        curr_tabs = []
        for t in tabs:
            nts = gen_next(t, partition, i + 1)
            curr_tabs.extend(nts)

        tabs = curr_tabs

    return tabs


def gen_next(curr, partition, el):
    next_parts = []

    for idx, row in enumerate(curr):
        if len(row) < partition[idx]: # there is space on this row
            if idx == 0 or (idx > 0 and len(row) < len(curr[idx - 1])):
                new_tab = copy.deepcopy(curr)
                new_tab[idx].append(el)
                next_parts.append(new_tab)

    if len(curr) < len(partition):
        new_tab = copy.deepcopy(curr)
        new_tab.append([el])
        next_parts.append(new_tab)

    return next_parts[::-1]

def pp(tableaux):
    for row in tableaux:
        fmt = '[{}]' * len(row)
        rowstr = fmt.format(*row)
        print(rowstr)

if __name__== '__main__':
    st = time.time()
    sty = gen_standard_tableaux((10,3,3))
    sty_time = time.time() - st

    st = time.time()
    yts = [YoungTableau(t, 30) for t in sty]
    yt_time = time.time() - st
    tot = sty_time + yt_time
    print(f'Elapsed: {tot:.2f} | Gen Shapes: {sty_time:.2f}s | Gen YT Obj: {yt_time:.2f}s')
    print(f'Prop obj: {yt_time/tot:.2f}')
    print(len(sty))
