def matcopy(matrix):
    num_rows = len(matrix)
    num_cols = len(matrix[0]) if matrix else 0
    new_matrix = [[0] * num_cols for _ in range(num_rows)]
    for i in range(num_rows):
        for j in range(num_cols):
            new_matrix[i][j] = matrix[i][j]
    
    return new_matrix

def matzeros(rows, cols):
    zero_matrix = [[0] * cols for _ in range(rows)]
    return zero_matrix

def veccopy(vector):
    copy = []
    for el in vector:
        copy.append(el)
    return copy

def vecsub(vector_1,vector_2):
    tmp = veccopy(vector_1)
    for i in range(len(tmp)):
        tmp[i]+=vector_2[i]
    return tmp

def vecones(length):
    vector = []
    for _ in range(length):
        vector.append(1)
    return vector

def veczeros(length):
    vector = []
    for _ in range(length):
        vector.append(0)
    return vector

def normalization(vector):
    i = 0
    for el in vector:
        i+=el**2
    return i**0.5

def diagToSquare(vector):
    size = len(vector)
    result_matrix = [[0] * size for _ in range(size)]  
    for i in range(size):
        result_matrix[i][i] = vector[i]  
    return result_matrix

