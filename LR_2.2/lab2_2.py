import math
import numpy as np
import matplotlib.pyplot as plt


def f_show1(x2):
  return math.cos(x2)/3

def f_show2(x2):
  return math.log(3*x2)

def d_phi1(x2):
    return -math.sin(x2)/3

def d_phi2(x1):
    return math.exp(x1)/3

def phi2(x1):
    return math.exp(x1)/3

def phi1(x2):
    return math.cos(x2)/3

def f1(x1, x2):
    return 3*x1 - math.cos(x2)

def f2(x1, x2):
    return 3*x2 - math.exp(x1)

def d_1fun_x1(x1, x2):
    return 3

def d_1fun_x2(x1, x2):
    return math.sin(x2)

def d_2fun_x1(x1, x2):
    return - math.exp(x1)

def d_2fun_x2(x1, x2):
    return 3

def det_a1(x1, x2):
    return f1(x1, x2) * d_2fun_x2(x1, x2) - f2(x1, x2) * d_1fun_x2(x1, x2)

def det_a2(x1, x2):
    return f2(x1, x2) * d_1fun_x1(x1, x2) - f1(x1, x2) * d_2fun_x1(x1, x2)

# Якобиан
def det_J(x1, x2):
    return d_1fun_x1(x1, x2) * d_2fun_x2(x1, x2) - d_2fun_x1(x1, x2) * d_1fun_x2(x1, x2)

def minus_vectors(a, b):
    n = len(a)
    res = np.zeros(n)

    for i in range(0, n):
        res[i] = a[i] - b[i]

    return res

def norma_vector(a):
    sum = 0
    n = len(a)
    for i in range(0, n):
        sum += a[i] * a[i]
        
    sum = math.sqrt(sum)
    return sum

def SimpleIterationsMethod(epsilon, n):
    xk = np.zeros(n)
    xk_tmp = np.zeros(n)

    xk[0] = 1.8
    xk[1] = 1
    x2a = 1.7
    x2b = 1.8
    x1a = 0.5
    x1b = 1
    #производные на участках мотононные, поэтому найдем максимум производных так:
    dphi1_max = max(d_phi1(x2a), d_phi1(x2b))
    dphi2_max = max(d_phi2(x1a), d_phi2(x1b))

    q = max(abs(dphi1_max), abs(dphi2_max))

    #проверка условия наличия решения
    if q > 1:
        print(f"Sufficient condition is not met : q = {q} > 1")
        return

    print(f"Sufficient condition is met : q = {q} < 1")

    counter = 0
    #Вычисляем корень по рекурентной формуле простых итераций
    while ((q / (1 - q)) * norma_vector(minus_vectors(xk, xk_tmp)) > epsilon):
        if counter != 0:
            xk[0] = xk_tmp[0]
            xk[1] = xk_tmp[1]

        counter = counter + 1
        xk_tmp[0] = phi1(xk[1])
        xk_tmp[1] = phi2(xk[0])

    print(f"x1 = {xk[0]}, x2 = {xk[1]}, on {counter} iteration")

def newtonMethod(x0, epsilon, n):
    
    xk = np.zeros(n)
    xk_tmp = np.zeros(n)
    detk = np.zeros(n)

    xk_tmp = x0

    counter = 0

    #Вычисляем корень по рекурентной формуле Ньютона
    while norma_vector(minus_vectors(xk, xk_tmp)) > epsilon:
        xk[0] = xk_tmp[0]
        xk[1] = xk_tmp[1]

        detk[0] = det_a1(xk[0], xk[1]) / det_J(xk[0], xk[1])
        detk[1] = det_a2(xk[0], xk[1]) / det_J(xk[0], xk[1])

        xk_tmp = minus_vectors(xk, detk)
        counter = counter + 1

    print(f"x1 = {xk[0]}, x2 = {xk[1]}, on {counter} iteration")

    
# МЕТОД ПРОСТЫХ ИТЕРАЦИЙ
    
# ГРАФИК
x = np.linspace(-2,2,100)
phi1_ = []
phi2_ = []
for x_ in x:
  phi1_.append(phi1(x_))
  phi2_.append(phi2(x_))
plt.plot(x,phi1_)
plt.plot(phi2_,x)
plt.xlabel('x2')
plt.ylabel('x1')
#положительное решение находится в квадрате 0 < x2 < 1, 0 < x1 < 0.5
#за начальное приближение возьмем точку (0.5,0.5)

epsilon = [0.1, 0.001, 0.0001, 0.00001]

# Посчитаем решения для различных точностей методом простых итераций
print("\n\t\tSimple iterations method\n")
for eps in epsilon:
    SimpleIterationsMethod(eps, 2)


# МЕТОД ГБЮТОНА

# ГРАФИК
x = np.linspace(0.05,1,100)
f1_ = []
f2_ = []
for x_ in x:
  f1_.append(f_show1(x_))
  f2_.append(f_show2(x_))
plt.plot(x,f1_)
plt.plot(x,f2_)
plt.xlabel('x2')
plt.ylabel('x1')

epsilon = [0.1, 0.001, 0.0001, 0.00001]

# Посчитаем решения для различных точностей методом Ньютона
print("\n\t\tNewton method\n")
for eps in epsilon:
    newtonMethod([0.5,0.5],eps, 2)
