#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Chenfei
@contact:chenfei.ye@foxmail.com
@version: 1.0.0
@file: surrogate_BOLD_create.py
@time: 2020/11/7 15:34
@function: create SC-informed surrogate BOLD signal
@ref: https://doi.org/10.1038/s41467-019-12765-7
"""

import numpy as np

def surrogate_BOLD_create(U, s, num_rand = 20):
    # Wsymm: input structural network
    # s: input BOLD signal
    # num_rand: number of surrogates
    # return: 2D array, #row = num_rand; #column = #s
    s_rand_ls = []
    for i in range(num_rand):
        np.random.seed(seed=i)
        PHIdiag = 2 * (np.random.randint(2, size= U.__len__()) - 0.5)
        PHI = np.diag(PHIdiag)
        XrandS = U @ PHI @ U.T @ s
        s_rand_ls.append(XrandS)
    return s_rand_ls
