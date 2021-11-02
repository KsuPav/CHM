"""lab1_3.ipynb
Лабораторная работа 1, задача 3, вариант 17
"""

import copy
from math import sqrt
import numpy as np


# МЕТОД ЗЕЙДЕЛЯ
def seidel(a, b, eps):
    n = len(a)
    beta = [0 for _ in range(n)]

    for i in range(n):
        beta[i] = b[i] / a[i][i]

    x = [beta[i] for i in range(n)]

    stop = False
    count = 0
    while not stop:
        x_new = copy.deepcopy(x)
        print("Iteration", count, "ans ", x)
        for i in range(n):
            tmp_1 = 0
            tmp_2 = 0
            for j in range(i):
                tmp_1 += a[i][j] * x_new[j] # Sum for x_(k+1)
            for j in range(i+1, n):
                tmp_2 += x[j] * a[i][j] # Sum for x_k

            x_new[i] = (b[i] - tmp_1 - tmp_2) / a[i][i]

        stop = sqrt(sum((x_new[i] - x[i]) ** 2 for i in range(n))) <= eps
        if not stop:
            x = x_new
        count = count + 1

    x = [round(x[i], 4) for i in range(n)]
    return x


# МЕТОД ПРОСТЫХ ИТЕРАЦИЙ
def simpleIteration(a, b, eps):
    n = len(b)
    beta = np.zeros(n)
    alpha = np.zeros((n, n))
    x = np.zeros(n)
    dx = np.zeros(n)
    # приводим систему к эквивалентному виду
    for i in range(n):
        beta[i] = b[i] / a[i][i]
        for j in range(n):
            if i != j:
                alpha[i][j] = - a[i][j] / a[i][i]
            else:
                alpha[i][j] = 0
        x[i] = beta[i]
        dx[i] = beta[i]

    # запускаем процесс, пока не выполнено условие остановки
    norm = 1
    norma_a = m_norm(n, alpha)
    count = 0
    while abs(norm) > eps:
        count = count + 1
        x = mv_mult(alpha, x)
        x = v_sum(x, beta)
        dx = v_minus(dx, x)
        norm = norma_a * v_norm(dx)/(1-norma_a)
        dx = x

    print(f"Решение найдено на итерации {count-1}:")
    print(x)

    
# СЛОЖЕНИЕ МАТРИЦ
def matrixsum(a, b):
    out = [[0] * len(a[0]) for _ in range(len(a))]
    for i in range(len(a)):
        for j in range(len(a[0])):
            out[i][j] = a[i][j] + b[i][j]
    return out


# УМНОЖЕНИЕ МАТРИЦЫ И ВЕКТОРА
def mv_mult(matrix, vector):
    size = len(vector)
    r = np.zeros(size)
    for i in range(size):
        for j in range(size):
            r[i] = r[i] + matrix[i][j] * vector[j]
    return r


# НОРМА МАТРИЦЫ (квадратный корень из суммы квадратов всех элементов матрицы)
def norma(x, y):
    return sqrt(sum((x[i] - y[i]) ** 2 for i in range(len(x))))


# НОРМА МАТРИЦЫ (максимальное из чисел, полученных при сложении всех элементов, взятых по модулю)
def m_norm(size, matrix):
    v = [abs(sum(matrix[i][j] for i in range(size))) for j in range(size)]
    return max(v)


# НОРМА ВЕКТОРА
def v_norm(vec):
    size = len(vec)
    return sum(abs(vec[i]) for i in range(size))


# СУММА ВЕКТОРОВ
def v_sum(x, y):
    return [x[i] + y[i] for i in range(len(x))]


# РАЗНОСТЬ ВЕКТОРОВ
def v_minus(x, y):
    return [x[i] - y[i] for i in range(len(x))]


# ВЫВОД МАТРИЦЫ
def show(a, n):
    for i in range(0, n):
        for j in range(0, n):
            print("\t", a[i][j], " ", end='')
        print("\n")


eps = 0.001
a = [[-19, 2, -1, -8],
    [2, 14, 0, -4],
    [6, -5, -20, -6],
    [-6, 4, -2, 15]]
b = [38, 20, 52, 43]
print("Метод Зейделя:")
seidel(a, b, eps)
print("\nМетод простых итераций:")
simpleIteration(a, b, eps)
