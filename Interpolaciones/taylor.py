import sympy as sp

class InterpolacionTaylor:
    def __init__(self, funcion, valor_x, num_derivadas, simbolo):
        self.funcion = funcion
        self.valor_x = valor_x
        self.num_derivadas = num_derivadas
        self.simbolo = simbolo
        self.polinomio = self._calcular_polinomio_taylor()

    def _calcular_polinomio_taylor(self):
        taylor = 0
        for i in range(self.num_derivadas + 1):
            derivada_i = sp.diff(self.funcion, self.simbolo, i)  # i-ésima derivada
            termino = (derivada_i.subs(self.simbolo, self.valor_x) / sp.factorial(i)) * (self.simbolo - self.valor_x)**i
            taylor += termino
        return sp.simplify(taylor)

    def graficar(self, rango=(-3, 5), puntos=100):
        funcion_original = sp.lambdify(self.simbolo, self.funcion, modules='sympy')
        funcion_taylor = sp.lambdify(self.simbolo, self.polinomio, modules='sympy')
        x_min, x_max = rango
        paso = (x_max - x_min) / (puntos - 1)
        valores_x_grafica = [x_min + i * paso for i in range(puntos)]
        valores_y_original = [funcion_original(x) for x in valores_x_grafica]
        valores_y_taylor = [funcion_taylor(x) for x in valores_x_grafica]
        lista = (valores_x_grafica,valores_y_original,valores_y_taylor)
        return lista

    def run(self):
        print("Polinomio de Taylor:",self.polinomio)
        self.graficar(rango=(-3, 5), puntos=200)
        lista = self.graficar()
        return lista