from math import sin
import time
import methods as meth
# Task A
def createMatrix(a1,a2,a3,N):
    A = []
    for i in range(N):
        row = []
        for j in range(N):
            if i == j:
                row.append(a1)
            elif abs(i - j) == 1:
                row.append(a2)
            elif abs(i - j) == 2:
                row.append(a3)
            else:
                row.append(0)
        A.append(row)
    return A

def createVector(N,f):
    b = []
    for n in range(N):
            b.append(sin((n+1) * (f + 1)))
    return b

# Task B
def jacobiMethod(A, b,iterations=3000):
    startTime = time.time()
    matA = meth.matcopy(A)
    vecB = meth.veccopy(b)
    k = 0 
    n = len(matA)
    x = meth.vecones(n)  # Initial guess for the solution vector
    residuals = []  # List to store residuals for each iteration
    while True:
        x_new = meth.vecones(n)  # Temporary storage for updated solution vector
        for i in range(n):
            sigma1 = sum(A[i][j] * x[j] for j in range(i))
            sigma2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (vecB[i] - sigma1-sigma2) / matA[i][i]
        
        # Calculate the residual vector
        residual = [sum(matA[i][j] * x_new[j] for j in range(n)) - vecB[i] for i in range(n)]
        residuals.append(meth.normalization(residual))
        
        # Check convergence
        if meth.normalization(residual)<1e-9:
            break
        
        # Update the solution vector
        x = x_new
        k += 1
    
    print("Jacobi's method")
    print('Time:', time.time() - startTime)
    print('Iterations:', k)
    print()
    return time.time() - startTime , residuals

        
def gaussSeidelMethod(A, b,iterations=3000):
    startTime = time.time()
    matA = meth.matcopy(A)
    vecB = meth.veccopy(b)
    k = 0
    n = len(matA)
    vecX = meth.vecones(n)
    residuals = []
    while True:
        x = vecX.copy()
        for i in range(n):
            sigma1 = sum(A[i][j] * vecX[j] for j in range(i))
            sigma2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            vecX[i] = (b[i] - sigma1 - sigma2) / A[i][i]
        
        residual = [sum(matA[i][j] * vecX[j] for j in range(n)) - vecB[i] for i in range(n)]
        residuals.append(meth.normalization(residual))
        if meth.normalization(residual) < 1e-9:
            break
        k += 1
    
    print("Gauss-Seidel's method")
    print("Time taken:", time.time() - startTime,"s")
    print("Iterations:", k)
    print()
    return time.time() - startTime, residuals

def luFactoryzationMethod(A, b):
    startTime = time.time()
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

    # Print results
    print("LU method")
    print('Time taken:', time.time() - startTime,'s')
    print("Residual norm:", meth.normalization(res))
    print()

    return time.time() - startTime

