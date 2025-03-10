import math
import matplotlib.pyplot as plt

class MetodoPotencia:
    def __init__(self, A, x0, tolerancia):
        self.A = A
        self.x0 = x0
        self.tolerancia = tolerancia
        self.autovalores, self.autovectores = self.metodo_potencia()

    def metodo_potencia(self, error = 1):
        if (len(self.A) == len(self.x0)):
            autovalores = []
            autovectores = []
            vector = self.x0
            while error > self.tolerancia:
                Ax = self.producto(self.A, vector) 
                norma = self.norma(Ax)
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
    
    def norma(self, vector):
        return math.sqrt(sum(x**2 for x in vector))
    
# Código principal:

matriz = [[4, -1, 1],
          [-1, 3, -2],
          [1, -2, 3]]

vector = [1, 0, 0]

tolerancia = 1e-8

mp = MetodoPotencia(matriz, vector, tolerancia)
print(f"Autovalor máximo: {mp.autovalores[-1]}")
print(f"Autovector normalizado: {mp.autovectores[-1]}")

# Comprobación:
print(f"\nAX = {mp.producto(mp.A, mp.autovectores[-1])}")
print(f"LambdaX = {[mp.autovalores[-1] * x for x in mp.autovectores[-1]]}")

plt.plot(mp.autovalores, marker='o')
plt.title('Convergencia del Autovalor Dominante')
plt.xlabel('Iteraciones')
plt.ylabel('Autovalor')
plt.grid()
plt.axhline(y=mp.autovalores[-1], color='r', linestyle='--', label='Valor final')
plt.legend()
plt.show()