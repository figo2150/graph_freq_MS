#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@version: 1.0.0
@file: surrogate_BOLD_create.py
@time: 2020/11/7 15:34
@functions: utilities for GSP_main.py
"""

import numpy as np
import pickle

def surrogate_BOLD_create(U, s, num_rand = 20):
    """
    function: create SC-informed surrogate BOLD signal
    Wsymm: input structural network
    s: input 2D or 1D BOLD signal
    num_rand: number of surrogates
    return: 3D or 2D array, #row = num_rand; #column = #s
    """
    if len(s.shape) == 1:
        s_rand = np.zeros(shape=(s.shape[0], num_rand))
    elif len(s.shape) == 2:
        s_rand = np.zeros(shape=(s.shape[0], s.shape[1], num_rand))
    else:
        ValueError('dimension error for functional signal')

    for i in range(num_rand):
        np.random.seed(seed=i)
        PHIdiag = 2 * (np.random.randint(2, size= U.__len__()) - 0.5)
        PHI = np.diag(PHIdiag)
        XrandS = U @ PHI @ U.T @ s
        if len(s.shape) == 1:
            s_rand[:,i] = XrandS
        else:
            s_rand[:,:,i] = XrandS
    return s_rand


def save_variable(v,filename):
    f=open(filename,'wb')
    pickle.dump(v,f)
    f.close()
    return filename
 
def load_variable(filename):
    f=open(filename,'rb')
    r=pickle.load(f)
    f.close()
    return r

# # calculate the median-split threshold 