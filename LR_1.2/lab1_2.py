# -*- coding: utf-8 -*-
"""lab1_2.ipynb

Лабораторная работа 1, задача 2, вариант 17
"""

import numpy as np

def runmethod(c, f, size):
    # Функция zeros() возвращает новый массив, заполненный нулями
    cn = np.zeros(size)
    fn = np.zeros(size)
    yn = np.zeros(size)
    xn = np.zeros(size)

    for i in range(0, size):
        if i == 0:
            yn[i] = c[i][i]
            cn[i] = -c[i][i + 1] / yn[i]
            fn[i] = f[i] / yn[i]
        elif i == size - 1:
            yn[i] = c[i][i] + c[i][i-1] * cn[i-1]
            fn[i] = (f[i] - c[i][i-1] * fn[i-1]) / yn[i]
        else:
            yn[i] = c[i][i] + c[i][i-1] * cn[i-1]
            cn[i] = -c[i][i+1] / yn[i]
            fn[i] = (f[i] - c[i][i-1] * fn[i-1]) / yn[i]

    xn[size-1] = fn[size-1]
    i = size-2
    
    while i >= 0:
        xn[i] = cn[i] * xn[i+1] + fn[i]
        i = i - 1

    return xn


a = [[-6, 5, 0, 0, 0],
      [-1, 13, 6, 0, 0],
      [0, -9, -15, -4, 0],
      [0, 0, -1, -7, 1],
      [0, 0, 0, 9, -18]]
b = [51, 100, -12, 47, -90]
print('Матрица А')
print(a)
print('Вектор b')
print(b)
print("Решение")
print(runmethod(a, b, len(b)))
