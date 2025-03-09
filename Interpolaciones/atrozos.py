import sympy as sp

class InterpolacionATrozos:
    def __init__(self, puntos):
        self.puntos = sorted(puntos)  # Ordenamiento de puntos.

    def interpolar_tramos_simbolico(self):
        x = sp.Symbol('x')
        tramos = []
        # Recorrer los segmentos y crear ecuaciones lineales para cada tramo.
        for i in range(len(self.puntos) - 1):
            x1, y1 = self.puntos[i]
            x2, y2 = self.puntos[i + 1]
            tramo = y1 + (y2 - y1) * (x - x1) / (x2 - x1)  # Ecuación de la recta.
            tramos.append((tramo, (x1, x2)))  # Guardar la ecuación y su intervalo.
        return tramos

    def interpolar(self, x_val):
        # Recorrer cada segmento.
        for i in range(len(self.puntos) - 1):
            x1, y1 = self.puntos[i]
            x2, y2 = self.puntos[i + 1]
            # Verificar si "x_val" está en el segmento actual.
            if x1 <= x_val <= x2:
                # Fórmula de la recta: y = y1 + (y2 - y1) * (x - x1) / (x2 - x1).
                return y1 + (y2 - y1) * (x_val - x1) / (x2 - x1)
        # Si no está en ningún segmento, devolver None (fuera del rango).
        return None

    def graficar(self, rango_x):
        x_grafica = []
        y_grafica = []
        for x in rango_x:
            y = self.interpolar(x)
            if y is not None:
                x_grafica.append(x)
                y_grafica.append(y)
        puntos_x = [p[0] for p in self.puntos]
        puntos_y = [p[1] for p in self.puntos]
        x_funcion = [x / 10.0 for x in range(int(rango_x[0] * 10), int(rango_x[-1] * 10) + 1)]
        y_funcion = [x**4 - 5*x**3 -2 for x in x_funcion]
        lista = (puntos_x,puntos_y,x_grafica,y_grafica,x_funcion,y_funcion)
        return lista
    
    def run(self):
        tramos = self.interpolar_tramos_simbolico()
        print("\nEcuaciones de los tramos de la Interpolación a Trozos:")
        for tramo, intervalo in tramos:
            print(f"Tramo en intervalo {intervalo}: {tramo}")
        rango = [x / 10.0 for x in range(-30, 51)]
        lista = self.graficar(rango)
        return lista