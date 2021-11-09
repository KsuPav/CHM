import matplotlib.pyplot as plt

def findInterval(xi, x0):
    for i in range(0,len(xi) - 1):
        if xi[i] <= x0 <= xi[i + 1]:
            return i

        
def firstDerivative(xi, yi, x0):
    i = findInterval(xi, x0)
    left = (yi[i + 1] - yi[i]) / (xi[i + 1] - xi[i])
    right = ((yi[i + 2] - yi[i + 1]) / (xi[i + 2] - xi[i + 1]) - left) / \
            (xi[i + 2] - xi[i]) * (2 * x0 - xi[i] - xi[i + 1])
    return left + right


def secondDerivative(xi, yi, x0):
    i = findInterval(xi, x0)
    left = (yi[i + 1] - yi[i]) / (xi[i + 1] - xi[i])
    right = 2 * ((yi[i + 2] - yi[i + 1]) / (xi[i + 2] - xi[i + 1]) - left) / \
            (xi[i + 2] - xi[i])
    return right


x0 = 0.2
xi = [-0.2, 0.0,	0.2,	0.4,	0.6]
yi = [-0.40136,	0.0,	0.40136,	0.81152,	1.2435]
print("Первая производная в точке = " + str(round(firstDerivative(xi, yi, x0), 6)))
print("Вторая производная в точке = " + str(round(secondDerivative(xi, yi, x0), 6)))


# Графики, чтобы убедиться в правильности результатов:
plt.figure(figsize=(10,8))
plt.scatter(xi,yi,linewidths=5)
plt.plot(xi,yi,c="yellow",linewidth=3)
plt.legend(["Заданная выборка","Попарное соединение двух соседних точек"])

"""Зависимость почти линейная, поэтому тангенс угла наклона касательной в 0.2 равен примерно двум, что совпадает с результатом. И так как функция почти линейна, то вторая производная примерно равна нулю, что видно из второго результата работы программы"""
