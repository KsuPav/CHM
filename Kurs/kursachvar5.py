import numpy as np
import matplotlib.pyplot as plt


def nonLinearFunction(x):
  return x**3/np.exp(x) + 3


def derrNonLinearFunction(x):
  return (3*x**2 - x**3)/np.exp(x)

x = np.linspace(-3,3,100)

plt.plot(x,nonLinearFunction(x))
plt.vlines(0,-3,5)
plt.hlines(0,-3,3)
plt.ylim(-3,5)
plt.xlim(-3,3)
plt.grid()
plt.legend(['Нелинейная функция'])


def homotopyMethod(x0):

  #число шагов для метода
  S = 10000

  #шаг
  h = 1/S

  #значение функции в начальной точке
  f0 = nonLinearFunction(x0)

  #первое значение параметра
  lk = 0

  xk = x0
  xkprev = x0
  
  #запускаем итеративный процесс, пока не параметр не станет равным 1
  for k in range(S):
    lk = lk + h
    dfk = derrNonLinearFunction(xk)
    xk = xkprev - h*f0/dfk
    xkprev = xk

  return xk

homotopyMethod(-1)

nonLinearFunction(-1.024887591013923)


def fun1(x1,x2):
  return x1*np.log(x2) - 1/x2


def dfun1x1(x1,x2):
  return np.log(x2)


def dfun1x2(x1,x2):
  return x1/x2 + 1/(x2**2)


def fun2(x1,x2):
  return x2 - np.exp(2*x1)


def dfun2x1(x1,x2):
  return -2*np.exp(2*x1)


def dfun2x2(x1,x2):
  return 1


def fun1_4show(x2):
  return (1/(x2*np.log(x2)))


def fun2_4show(x2):
  return np.log(x2)/2


def jacobianDet(x1,x2):
  return dfun1x1(x1,x2)*dfun2x2(x1,x2) - dfun1x2(x1,x2)*dfun2x1(x1,x2)


def jacobianObr(x1,x2):
  return [[dfun2x2(x1,x2)/jacobianDet(x1,x2),-dfun1x2(x1,x2)/jacobianDet(x1,x2)],
          [-dfun2x1(x1,x2)/jacobianDet(x1,x2),dfun1x1(x2,x2)/jacobianDet(x1,x2)]]


x = np.linspace(0.25,4,100)

plt.plot(x,fun1_4show(x))
plt.plot(x,fun2_4show(x))
plt.xlim(1,4)
plt.ylim(0.25,1)
plt.xlabel('x2')
plt.ylabel('x1')

#начальное приближение - (0.45,2.5)

print(fun1(0.44,2.46))
print(fun2(0.44,2.46))

def matrixMultVector(m,v):
  retVal = [0,0]
  for m_r in range(0,len(m)):
    for m_c in range(0,len(m)):
      retVal[m_r] += m[m_r][m_c]*v[m_c] 
  return retVal

def homotopyMethodSystem(x0):

  #число шагов для метода
  S = 1000

  #шаг
  h = 1/S

  #значение функции в начальной точке
  f0 = [fun1(x0[0],x0[1]),fun2(x0[0],x0[1])]
  print(f0)

  #первое значение параметра
  lk = 0

  xk = x0
  xkprev = x0
  
  #запускаем итеративный процесс, пока не параметр не станет равным 1
  for k in range(S):
    lk = lk + h
    detJ = jacobianDet(xk[0],xk[1])
    ftemp = [(h)*f0[0],(h)*f0[1]]
    ftemp = matrixMultVector(jacobianObr(xk[0],xk[1]),ftemp)
    xk = [xkprev[0] - ftemp[0],xkprev[1] - ftemp[1]]
    xkprev = xk

  return xk

solve = homotopyMethodSystem([0.45,2.5])

print(fun1(solve[0], solve[1]))
print(fun2(solve[0], solve[1]))
