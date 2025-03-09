import sympy as sp

class InterpolacionLagrange:
    def __init__(self, puntos):
        self.puntos = puntos
        self.x_simbolo = sp.Symbol('x')  # Variable simbólica para sympy

    def calcular_polinomio_sympy(self):
        polinomio = 0
        for i in range(len(self.puntos)):
            termino = self.puntos[i][1]  # Valor y de i
            for j in range(len(self.puntos)):
                if i != j:
                    termino *= (self.x_simbolo - self.puntos[j][0]) / (self.puntos[i][0] - self.puntos[j][0])
            polinomio += termino
        return str(sp.simplify(polinomio))

    def calcular_polinomio(self, x_val):
        resultado = 0
        for i in range(len(self.puntos)):
            termino = self.puntos[i][1]  # Valor y de i
            for j in range(len(self.puntos)):
                if i != j:
                    termino *= (x_val - self.puntos[j][0]) / (self.puntos[i][0] - self.puntos[j][0])
            resultado += termino
        return resultado

    def graficar(self, rango_x):
        x_grafica = []
        y_grafica = []
        # Evaluar el polinomio en cada valor de rango_x
        for x in rango_x:
            x_grafica.append(x)
            y_grafica.append(self.calcular_polinomio(x))
        # Graficar los puntos originales
        puntos_x = [p[0] for p in self.puntos]
        puntos_y = [p[1] for p in self.puntos]
        lista = (puntos_x,puntos_y,x_grafica,y_grafica)
        return lista

    def run(self):
        print("\nPolinomio de Lagrange:", self.calcular_polinomio_sympy())
        rango = [x / 10.0 for x in range(-25, 50)]
        lista = self.graficar(rango)
        return lista