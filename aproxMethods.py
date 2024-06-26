import os
import csv
import methods as meth
import equationSolving as initMeth
import matplotlib.pyplot as plt
import numpy as np

def chebyshev_nodes(a, b, n):
    return [0.5 * (a + b) + 0.5 * (b - a) * np.cos(i* np.pi  / (n-1) ) for i in range(n)]

def LagrangeMethod(k, use_chebyshev=False):
    def interpolation_function(points):
        def f(x):
            result = 0
            n = len(points)
            for i in range(n):
                xi, yi = points[i]
                base = 1
                for j in range(n):
                    if i != j:
                        xj, yj = points[j]
                        base *= (float(x) - float(xj)) / (float(xi) - float(xj))
                result += float(yi) * base
            return result
        return f

    for file in os.listdir('./data'):
        with open(f'./data/{file}', 'r') as f:
            data = list(csv.reader(f))
        
        data = [(float(x), float(y)) for x, y in data[1:]]
        if use_chebyshev:
            a, b = data[0][0], data[-1][0]
            nodes = chebyshev_nodes(a, b, len(data) // k)
            interpolation_data = [(node, next(y for x, y in data if x >= node)) for node in nodes]
        else:
            interpolation_data = data[1::k]
        
        F = interpolation_function(interpolation_data)

        distance, height, interpolated_height = [], [], []
        for x, y in data:
            distance.append(float(x))
            height.append(float(y))
            interpolated_height.append(F(float(x)))

        train_distance, train_height = [], []
        for x, y in interpolation_data:
            train_distance.append(float(x))
            train_height.append(F(float(x)))

        plt.semilogy(distance, height, 'r.', label='pełne dane')
        plt.semilogy(distance, interpolated_height, color='blue', label='funkcja interpolująca')
        plt.semilogy(train_distance, train_height, 'g.', label='dane do interpolacji')
        plt.legend()
        plt.ylabel('Wysokość')
        plt.xlabel('Odległość')
        plt.title('Przybliżenie interpolacją Lagrange\'a, ' + str(len(interpolation_data)) + ' punkty(ów)')
        plt.suptitle(file)
        plt.grid()
        plt.show()

def splineMethod(k):
    def interpolation_function(points):
        def calculate_params():
            n = len(points)
            A = meth.matzeros(4 * (n - 1), 4 * (n - 1))
            b = meth.veczeros(4 * (n - 1))

            for i in range(n - 1):
                x, y = points[i]
                row = meth.veczeros(4 * (n - 1))
                row[4 * i + 3] = 1
                A[4 * i + 3] = row
                b[4 * i + 3] = float(y)

            for i in range(n - 1):
                x1, y1 = points[i + 1]
                x0, y0 = points[i]
                h = float(x1) - float(x0)
                row = meth.veczeros(4 * (n - 1))
                row[4 * i] = h ** 3
                row[4 * i + 1] = h ** 2
                row[4 * i + 2] = h
                row[4 * i + 3] = 1
                A[4 * i + 2] = row
                b[4 * i + 2] = float(y1)

            for i in range(n - 2):
                x1, y1 = points[i + 1]
                x0, y0 = points[i]
                h = float(x1) - float(x0)
                row = meth.veczeros(4 * (n - 1))
                row[4 * i] = 3 * h ** 2
                row[4 * i + 1] = 2 * h
                row[4 * i + 2] = 1
                row[4 * (i + 1) + 2] = -1
                A[4 * i] = row
                b[4 * i] = 0.0

            for i in range(n - 2):
                x1, y1 = points[i + 1]
                x0, y0 = points[i]
                h = float(x1) - float(x0)
                row = meth.veczeros(4 * (n - 1))
                row[4 * i] = 6 * h
                row[4 * i + 1] = 2
                row[4 * (i + 1) + 1] = -2
                A[4 * (i + 1) + 1] = row
                b[4 * (i + 1) + 1] = 0.0

            row = meth.veczeros(4 * (n - 1))
            row[1] = 2
            A[1] = row
            b[1] = 0.0

            row = meth.veczeros(4 * (n - 1))
            x1, y1 = points[-1]
            x0, y0 = points[-2]
            h = float(x1) - float(x0)
            row[1] = 2
            row[-4] = 6 * h
            A[-4] = row
            b[-4] = 0.0

            result = initMeth.luFactoryzationMethod(A, b)
            return result

        params = calculate_params()

        def f(x):
            param_array = []
            row = []
            for param in params:
                row.append(param)
                if len(row) == 4:
                    param_array.append(row.copy())
                    row.clear()

            for i in range(1, len(points)):
                xi, yi = points[i - 1]
                xj, yj = points[i]
                if float(xi) <= x <= float(xj):
                    a, b, c, d = param_array[i - 1]
                    h = x - float(xi)
                    return a * h ** 3 + b * h ** 2 + c * h + d

            if x < float(points[0][0]):
                return param_array[0][3]
            elif x > float(points[-1][0]):
                return param_array[-1][3]

            return None

        return f

    for file in os.listdir('./data'):
        with open(f'./data/{file}', 'r') as f:
            data = list(csv.reader(f))

        data = data[1:]
        shift = (-1) * (len(data) % k)
        if shift != 0:
            interpolation_data = data[:shift:k]
        else:
            interpolation_data = data[::k]

        if interpolation_data[-1] != data[-1]:
            interpolation_data.append(data[-1])

        interpolation_data = [(float(x), float(y)) for x, y in interpolation_data]
        F = interpolation_function(interpolation_data)

        distance = []
        height = []
        interpolated_height = []
        for point in data:
            x, y = point
            distance.append(float(x))
            height.append(float(y))
            interpolated_height.append(F(float(x)))

        train_distance = []
        train_height = []
        for point in interpolation_data:
            x, y = point
            train_distance.append(float(x))
            train_height.append(float(y))

        plt.plot(distance, height, 'r.', label='pełne dane')
        plt.plot(distance, interpolated_height, color='blue', label='funkcja interpolująca')
        plt.plot(train_distance, train_height, 'g.', label='dane do interpolacji')
        plt.legend()
        plt.ylabel('Wysokość')
        plt.xlabel('Odległość')
        plt.title('Przybliżenie interpolacją Splajnami, ' + str(len(interpolation_data)) + ' punkty(ów)')
        plt.suptitle(file)
        plt.grid()
        plt.show()