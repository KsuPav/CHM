import numpy as np
import math

"""Лабораторная работа 1, задача 4, вариант 17"""


# МЕТД ВРАЩЕНИЙ
def rotation(a, eps=0.01):
    n = len(a)
    ak = [row.copy() for row in a]

    u = [[0. if i != j else 1. for i in range(n)] for j in range(n)]

    cov = False

    while not cov:
        ik, jk = 0, 1

        # Выбирается максимальный по модулю недиагональный элемент
        for i in range(n - 1):
            for j in range(i + 1, n):
                if abs(ak[i][j]) > abs(ak[ik][jk]):
                    ik, jk = i, j
        if ak[ik][ik] == ak[jk][jk]:
            phi = math.pi / 4
        else:
            phi = math.atan(2 * a[ik][jk] / (a[ik][ik] - a[jk][jk])) * 0.5

        uk = [[0. if i != j else 1. for i in range(n)] for j in range(n)]

        uk[ik][jk] = math.sin(phi)
        uk[jk][ik] = -math.sin(phi)

        uk[ik][ik] = math.cos(phi)
        uk[jk][jk] = math.cos(phi)

        # домножаем матрицу а слева на u^T и справа на u^T
        tmp = multi(uk, ak)
        uk[ik][jk], uk[jk][ik] = uk[jk][ik], uk[ik][jk]

        ak = multi(tmp, uk)
        u = multi(u, uk)

        count = 0

        for i in range(n - 1):
            for j in range(i + 1, n):
                count += ak[i][j] ** 2

        average = math.sqrt(count)
        if average < eps:
            cov = True

    return [ak[i][i] for i in range(n)], u


# УМНОЖЕНИЕ МАТРИЦ
def multi(m1, m2):
    sum = 0  # сумма
    tmp = []  # временная матрица
    ans = []  # конечная матрица
    row1 = len(m1)  # количество строк в первой матрице
    col1 = len(m1[0])  # Количество столбцов в 1
    row2 = col1  # и строк во 2ой матрице
    col2 = len(m2[0])  # количество столбцов во 2ой матрице
    for k in range(0, row1):
        for j in range(0, col2):
            for i in range(0, col1):
                sum = sum + m1[k][i] * m2[i][j]
            tmp.append(sum)
            sum = 0
        ans.append(tmp)
        tmp = []
    return ans


# СЛОЖЕНИЕ МАТРИЦ
def matrixsum(a, b):
    out = [[0] * len(a[0]) for _ in range(len(a))]
    for i in range(len(a)):
        for j in range(len(a[0])):
            out[i][j] = a[i][j] + b[i][j]
    return out


# ВЫВОД МАТРИЦЫ
def show(a, n):
    for i in range(0, n):
        for j in range(0, n):
            print("\t", a[i][j], " ", end='')
        print("\n")


eps = 0.01
n = 3

a = [[5, -3, -4], [-3, -3, 4], [-4, 4, 0]]

print("Симметрическая матрица:")
show(a, n)

x, u = rotation(a, eps)
print('x:\n', x)
print('u:\n')
show(u, len(u))

print("Проверка с помощью linalg:")
x, u = np.linalg.eig(a)
print('x:\n', x)
print('u:\n')
show(u, len(u))
