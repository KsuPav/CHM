def getCofactor(a, tmp, p, q, n):
  i = 0
  j = 0
  #копируем в новую матрицу все элементы исходной без строки и столбца
  for row in range(n):
    for col in range(n):
      if row != p and col != q:
        tmp[i][j] = a[row][col]
        j += 1
        #переход к первому столбцу и на следующую строку
        if j == n - 1:
          j = 0
          i += 1

# ДЕТЕРМИНАНТ
def determinant(a, n):
    d = 0
    if n == 1:
        return a[0][0]

    tmp = [[0] * n for _ in range(n)]  # Cofactors
    sign = 1  # коэффициент при слагаемом (-1)^k

    for i in range(n):
        getCofactor(a, tmp, 0, i, n)
        d = d + sign * a[0][i] * determinant(tmp, n - 1)

        sign = -sign

    return d


# СОПРЯЖЕННАЯ МАТРИЦА
def adjoin(a, n):
    adj = [[0] * n for _ in range(n)]
    if n == 1:
        adj[0][0] = 1
        return

    tmp = [[0] * n for _ in range(n)]  # Cofactors
    for i in range(n):
        for j in range(n):
            getCofactor(a, tmp, i, j, n)  # Cofactor a[i][j]

            # sign = (-1)**(i+j)
            if (i + j) % 2 == 0:
                sign = 1
            else:
                sign = -1

            adj[j][i] = sign * (determinant(tmp, n - 1))  # (-1^(i+j)*|алгебр дополнение|)

    return adj

  
# ОБРАТНАЯ МАТРИЦА
def inverse(a, b, n):
    det = determinant(a, n)
    if det == 0:
        print("Матрица вырождена")
        return False

    adj = adjoin(a, n)

    for i in range(n):
        for j in range(n):
            b[i][j] = adj[i][j] / det

    return True

  
# ТРАНСПОНИРОВАНИЕ МАТРИЦЫ
def transpose(a, n):
    b = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            b[i][j] = a[j][i]
    return b

# ПЕРЕМНОЖЕНИЕ МАТРИЦ
def multi(M1, M2):
    sum = 0
    tmp = []
    ans = []

    row1 = len(M1) 
    col1 = len(M1[0]) 
    row2 = col1
    col2 = len(M2[0])
    for k in range(0, row1):
        for j in range(0, col2):
            for i in range(0, col1):
                sum = round(sum + M1[k][i] * M2[i][j], 8)
            tmp.append(sum)
            sum = 0
        ans.append(tmp)
        tmp = []
    return ans

  
# 
def lup_solve(l, u, pi, b, n):
    x = [0 for i in range(n)]
    y = [0 for i in range(n)]

    for i in range(n):
        summ = 0
        for j in range(i):
            summ += l[i][j] * y[j]

        y[i] = b[pi[i]] - summ

    for i in range(n - 1, -1, -1):
        sec_summ = 0
        for j in range(i + 1, n):
            sec_summ += u[i][j] * x[j]

        x[i] = (y[i] - sec_summ) / u[i][i]

    x = [round(x[i], 5) for i in range(len(x))]
    return x


# 
def lupdecompose(a, n):
    pi = [i for i in range(n)]

    for k in range(n):
        p = 0
        for i in range(k, n):
            if abs(a[i][k]) > p:
                p = abs(a[i][k])
                tmp_k = i

        pi[k], pi[tmp_k] = pi[tmp_k], pi[k]

        for i in range(n):
            a[k][i], a[tmp_k][i] = a[tmp_k][i], a[k][i]
        for i in range(k + 1, n):
            a[i][k] = a[i][k] / a[k][k]
            for j in range(k + 1, n):
                a[i][j] = a[i][j] - a[i][k] * a[k][j]
    return pi    
    
    
# 
def get_lu(a):
    n = len(a)
    l = [[0] * n for i in range(0, n)]
    u = [[0] * n for i in range(0, n)]

    for i in range(n):
        l[i][i] = 1
        for j in range(n):
            if j < i:
                l[i][j] = a[i][j]
            else:
                u[i][j] = a[i][j]
    return l, u

  
# 
def roundMatrix(a, after):
  retVal = [[0] * len(a) for _ in range(len(a[0]))]

  for i in range(0,len(a)):
    for j in range(0,len(a[0])):
      retVal[i][j] = round(a[i][j],after)

  return retVal


# ВЫВОД МАТРИЦЫ
def show(a, n):
    for i in range(0, n):
        for j in range(0, n):
            print(" ", a[i][j], " ", end="")
        print("\n")
