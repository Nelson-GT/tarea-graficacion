# Optimización del Diseño de un Puente mediante Análisis Estructural.   

import math
import matplotlib.pyplot as plt

class Determinante:
    def __init__(self):
        pass
    
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

class Inversa:
    def __init__(self):
        pass

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
    
class NumeroCondicion(Determinante, Inversa):
    def __init__(self, A):
        self.A = A
        Determinante.__init__(self)
        self.detA = self.determinante(A)
        Inversa.__init__(self)
        self.condA = self.calcular_condA()

    def calcular_condA(self):
        if self.detA == 0:
            print("La matriz A no es invertible.")
            return None
        else:
            invA = self.inversa(self.A)
            return self.norma_matriz(self.A) * self.norma_matriz(invA)
    
    def copiar_matriz(self, matriz):
        return [fila[:] for fila in matriz]
    
    def norma_matriz(self, matriz):
        return math.sqrt(sum(x**2 for fila in matriz for x in fila))
    
class MetodoPotencia:
    def __init__(self, A, x0, tolerancia):
        self.A = A
        self.x0 = x0
        self.tolerancia = tolerancia
        self.autovaloresmax, self.autovectoresmax = self.metodo_potencia()

    def metodo_potencia(self, error = 1):
        if (len(self.A) == len(self.x0)):
            autovalores = []
            autovectores = []
            vector = self.x0
            while error > self.tolerancia:
                Ax = self.producto(self.A, vector) 
                norma = self.norma_vectormp(Ax)
                autovalores.append(norma)
                vector = [x / norma for x in Ax]
                autovectores.append(vector)
                if (len(autovalores) > 1):
                    error = abs(autovalores[-1] - autovalores[-2])
            return autovalores, autovectores
        else:
            print("Error de dimensiones")
        
    def producto(self, matriz, vector):
        return [sum(a * b for a, b in zip(fila, vector)) for fila in matriz]
    
    def norma_vectormp(self, vector):
        return math.sqrt(sum(x**2 for x in vector))
    
class MetodoPotenciaInversa(Determinante, Inversa):
    def __init__(self, A, x0, tolerancia):
        self.A = A
        self.x0 = x0
        self.tolerancia = tolerancia
        Determinante.__init__(self)
        self.detA = self.determinante(A)
        if self.detA != 0:
            Inversa.__init__(self)
            self.invA = self.inversa(A)
            self.autovaloresmin, self.autovectoresmin = self.metodo_potenciainv()
        else:
            print("La matriz no es invertible")
    
    def copiar_matriz(self, matriz):
        return [fila[:] for fila in matriz]

    def metodo_potenciainv(self, error = 1):
        if (len(self.A) == len(self.x0)):
            autovalores = []
            autovectores = []
            vector = self.x0
            while error > self.tolerancia:
                invAx = self.producto(self.invA, vector) 
                norma = self.norma_vectormpi(invAx)
                autovalores.append(1/norma)
                vector = [x / norma for x in invAx]
                autovectores.append(vector)
                if (len(autovalores) > 1):
                    error = abs(autovalores[-1] - autovalores[-2])
            return autovalores, autovectores
        else:
            print("Error de dimensiones")
        
    def producto(self, matriz, vector):
        return [sum(a * b for a, b in zip(fila, vector)) for fila in matriz]
    
    def norma_vectormpi(self, vector):
        return math.sqrt(sum(x**2 for x in vector))
    
class MatrizRigidez(NumeroCondicion, MetodoPotencia, MetodoPotenciaInversa):
    def __init__(self, A):
        NumeroCondicion.__init__(self, A)
        MetodoPotencia.__init__(self, A, [1, 0, 0], 1e-12)
        MetodoPotenciaInversa.__init__(self, A, [1, 0, 0], 1e-12)

# Código Principal:

# Matriz de rigidez: Describe la relación entre fuerzas aplicadas y desplazamiento de puntos clave del puente.
matrizrigidez = [[4, 8, 9], 
                [8, 8, 5], 
                [9, 5, 1]]

mc = MatrizRigidez(matrizrigidez)

# Número de condición de la matriz de rigidez.
print("--- Resistencia de la estructura a fuerzas aplicadas ---")
print(f"Número de condición: {mc.condA}")

# Autovalor dominante de la matriz de rigidez.
print("\n--- Autovalor dominante: Puntos de mayor rigidez de la estructura ---")
print(f"\nAutovalor máximo: {mc.autovaloresmax[-1]}")
print(f"Autovector normalizado (dirección de mayor resistencia estructural): {mc.autovectoresmax[-1]}")

# Comprobación:
print(f"\nComprobación\nAX = {mc.producto(mc.A, mc.autovectoresmax[-1])}")
print(f"LambdaX = {[mc.autovaloresmax[-1] * x for x in mc.autovectoresmax[-1]]}")

plt.plot(mc.autovaloresmax, marker='o')
plt.title('Convergencia del Autovalor Dominante: Puntos de mayor rigidez estructural')
plt.xlabel('Iteraciones')
plt.ylabel('Autovalor')
plt.grid()
plt.axhline(y=mc.autovaloresmax[-1], color='r', linestyle='--', label='Valor final')
plt.legend()
plt.show()

# Autovalor mínimo de la matriz de rigidez.
if (mc.detA != 0):
    print("\n--- Autovalor mínimo: Puntos de menor rigidez de la estructura ---")
    print(f"\nAutovalor mínimo: {mc.autovaloresmin[-1]}")
    print(f"Autovector normalizado (dirección de mayor vulnerabilidad estructural): {mc.autovectoresmin[-1]}")

    # Comprobación:
    print(f"\nComprobación\ninvAX = {mc.producto(mc.invA, mc.autovectoresmin[-1])}")
    print(f"LambdaX = {[(1/mc.autovaloresmin[-1]) * x for x in mc.autovectoresmin[-1]]}")

    plt.plot(mc.autovaloresmin, marker='o')
    plt.title('Convergencia del Autovalor Mínima: Puntos de menor rigidez estructural')
    plt.xlabel('Iteraciones')
    plt.ylabel('Autovalor')
    plt.grid()
    plt.axhline(y=mc.autovaloresmin[-1], color='r', linestyle='--', label='Valor final')
    plt.legend()
    plt.show()