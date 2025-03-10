import math
import random
import matplotlib.pyplot as plt

class NumeroCondicion:
    def __init__(self, A):
        self.A = A
        self.condA = self.calcular_condA()

    def calcular_condA(self):
        detA = self.determinante(self.A)
        if detA == 0:
            print("La matriz A no es invertible.")
            return None
        else:
            invA = self.inversa(self.A)
            return self.norma(self.A) * self.norma(invA)
    
    def determinante(self, matriz):
        if len(matriz) == 2:
            return ((matriz[0][0] * matriz[1][1]) - (matriz[0][1] * matriz[1][0]))
        determinante = 0
        for c in range(len(matriz)):
            # Crear una submatriz. Expansión por cofactores a partir de la fila 1.
            submatriz = [fila[:c] + fila[c+1:] for fila in matriz[1:]]
            # Calcular el determinante.
            determinante += ((-1) ** c) * matriz[0][c] * self.determinante(submatriz)
        return determinante

    def inversa(self, matriz):
        n = len(matriz)
        A = self.copiar_matriz(matriz)
        identidad = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
        for i in range(n):
            A[i] += identidad[i]
        for i in range(n):
            if abs(A[i][i]) < 1e-12:
                for j in range(i + 1, n):
                    if abs(A[j][i]) > abs(A[i][i]):
                        A[i], A[j] = A[j], A[i]
                        break
            pivote = A[i][i]
            A[i] = [x / pivote for x in A[i]]
            for j in range(n):
                if j != i:
                    factor = A[j][i]
                    A[j] = [a - factor * b for a, b in zip(A[j], A[i])]
        return [fila[n:] for fila in A]
    
    def copiar_matriz(self, matriz):
        return [fila[:] for fila in matriz]
    
    def norma(self, matriz):
        return math.sqrt(sum(x**2 for fila in matriz for x in fila))
    
# Código principal:
x = []
y = []
n = 2

for i in range(5):
    nc = NumeroCondicion([[random.uniform(-100, 100) for _ in range(n)] for _ in range(n)])
    if nc.condA != None:
        x.append(n)
        y.append(nc.condA)
    n += 1
    print(f"Número de condición de la matriz {i + 1}: {nc.condA}")

plt.plot(x, y, marker = 'o')
plt.xlabel('Tamaño de la matriz (n x n)')
plt.ylabel('Número de condición de la matriz')
plt.title('cond(A) vs n')
plt.grid()
plt.show()