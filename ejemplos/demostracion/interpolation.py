# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 21:20:15 2016

@author: cosmoscalibur
"""

def lagrange_polynomial(x_points, x_index, x_eval):
    N = len(x_points)
    cont = 0
    poly = 1
    while cont < N:
        if cont != x_index:
            poly *= (x_eval - x_points[cont]) / (x_points[x_index] - x_points[cont])
        cont += 1
    return poly

def lagrange_interpolation(x_points, y_points, x_eval):
    N = len(x_points)
    cont = 0
    poly = 0
    while cont < N:
        if hasattr(y_points, '__iter__'):
            poly += y_points[cont] * lagrange_polynomial(x_points, cont, x_eval)
        else:
            poly += y_points(x_points[cont]) * lagrange_polynomial(x_points, cont, x_eval)
        cont += 1
    return poly