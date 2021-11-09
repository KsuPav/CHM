import numpy as np


def func(x):
    return 1/(256-x**4)


def getX(x0, x, step):
    size = 1 + int((x - x0)/step)
    return np.linspace(x0, x, size)


def getY(x):
    return [func(i) for i in x]


def rectangle(x, h):
    return h * sum([func((x[i] + x[i + 1]) / 2) for i in range(len(x) - 1)])


def trapeze(x, h):
    y = getY(x)
    return h * (y[0] / 2 + sum([y[i] for i in range(1, len(y) - 2)]) + y[len(y) - 1] / 2)


def simpson(x, h):
    y = getY(x)
    return h / 3 * (y[0] +
                    sum([4 * y[i] for i in range(1, len(y) - 1, 2)]) +
                    sum([2 * y[i] for i in range(2, len(y) - 2, 2)]) +
                    y[len(y) - 1])


def rungerombergerror(res):
    f1 = res[0][0]
    h1 = res[0][1]
    f2 = res[1][0]
    h2 = res[1][1]
    k = max(h1/h2, h2/h1)
    return abs(f1-f2)/(k**2 - 1)


x0 = 0
xn = 2
h = [1.0, 0.5]
true_value = (np.log(3) + 2*np.arctan(0.5))/256


def printRes(method, methodName, h, x0, xn):
    print('\tМетод ' + methodName + ':')
    res = []
    for h_i in h:    
        print('Шаг h =' + str(h_i))
        x = getX(x0, xn, h_i)
        y = getY(x)

        
        res_rec = method(x, h_i)
        res.append([res_rec, h_i])
        print("Значение интеграла = " + str(res_rec))
    return res, res_rec


def getErrors(res, res_rec, true_value):
  print("\nОшибка Рунге-Ромберга " + str(rungerombergerror(res)))
  print("Ошибка абсолютная = " + str(abs(res_rec - true_value)))
  print('\n\n')


res, res_rec = printRes(rectangle,"rectangle",h,x0,xn)
getErrors(res, res_rec, true_value)

res, res_rec = printRes(trapeze,"trapeze",h,x0,xn)
getErrors(res, res_rec, true_value)

res, res_rec = printRes(simpson,"simpson",h,x0,xn)
getErrors(res, res_rec, true_value)
