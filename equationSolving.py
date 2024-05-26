import methods as meth

def luFactoryzationMethod(A, b):
    
    n = len(A)

    matA = meth.matcopy(A)
    mat_l = meth.diagToSquare(meth.vecones(n))
    mat_u = meth.matzeros(n, n)

    vecB = meth.veccopy(b)
    vecX = meth.vecones(n)
    vecY = meth.veczeros(n)

    # LUx = b
    for j in range(n):
        # Calculate U elements
        for i in range(j + 1):
            mat_u[i][j] += matA[i][j]
            for k in range(i):
                mat_u[i][j] -= mat_l[i][k] * mat_u[k][j]
        
        # Calculate L elements
        for i in range(j + 1, n):
            for k in range(j):
                mat_l[i][j] -= mat_l[i][k] * mat_u[k][j]
            mat_l[i][j] += matA[i][j]
            mat_l[i][j] /= mat_u[j][j]

    # Ly = b
    for i in range(n):
        value = vecB[i]
        for j in range(i):
            value -= mat_l[i][j] * vecY[j]
        vecY[i] = value / mat_l[i][i]

    # Ux = y
    for i in range(n - 1, -1, -1):
        value = vecY[i]
        for j in range(i + 1, n):
            value -= mat_u[i][j] * vecX[j]
        vecX[i] = value / mat_u[i][i]

    # Calculate residual
    res = [sum(matA[i][j] * vecX[j] for j in range(n)) - vecB[i] for i in range(n)]

    return vecX

