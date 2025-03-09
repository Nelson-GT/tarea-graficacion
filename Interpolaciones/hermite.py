import sympy as sp

class InterpolacionHermite:
    def __init__(self, valores_x, valores_y, derivadas):
        self.valores_x = valores_x
        self.valores_y = valores_y
        self.derivadas = derivadas
        self.polinomio = self._calcular_polinomio_hermite()

    def _calcular_polinomio_hermite(self):
        x = sp.symbols('x')
        n = len(self.valores_x)
        H = 0
        for i in range(n):
            # Calcula el producto base LaGrange.
            h_i = 1
            for j in range(n):
                if i != j:
                    h_i *= (x - self.valores_x[j]) / (self.valores_x[i] - self.valores_x[j])
            derivada_h_i = sp.diff(h_i, x)
            # Construye el polinomio de Hermite.
            H += (1 - 2 * derivada_h_i.subs(x, self.valores_x[i]) * (x - self.valores_x[i])) * h_i**2 * self.valores_y[i]
            H += (x - self.valores_x[i]) * h_i**2 * self.derivadas[i]
        return sp.simplify(H)

    def graficar(self, puntos=100):
        x = sp.symbols('x')
        funcion_hermite = sp.lambdify(x, self.polinomio, modules='sympy')
        # Rango de la funcion: Xmenor - 1, Xmayor + 3 (-2-1 = -3; 2+3 = 5).
        x_min, x_max = min(self.valores_x) - 1, max(self.valores_x) + 3
        paso = (x_max - x_min) / (puntos - 1)
        valores_x_grafica = [x_min + i * paso for i in range(puntos)]
        valores_y_grafica = [funcion_hermite(x_val) for x_val in valores_x_grafica]
        lista = (valores_x_grafica,valores_y_grafica,self.valores_x,self.valores_y)
        return lista
    
    def run(self):
        print("\nPolinomio de Hermite: ",self.polinomio)
        lista = self.graficar()
        return lista