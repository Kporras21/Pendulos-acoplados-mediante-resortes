import numpy as np
import sympy as sp
from scipy.constants import g

def normal_frec(k1, k2, m, L):
    # Definimos la variable simbólica
    w2 = sp.symbols('w2')

    # Definimos la matriz simbólica B con la variable w2
    B = sp.Matrix([
        [m * g * L + (L**2) * (k1 + k2) - w2 * m * (L**2), -k2 * (L**2)],
        [-k2 * (L**2), m * g * L + (L**2) * k2 - w2 * m * (L**2)]
    ])

    # Calculamos el determinante de la matriz B
    det_B = B.det()
    print(f'Su matriz es:\n{B}\n\nEl determinante de esta matriz es: \n{det_B}')

    # Resolvemos la ecuación det(B) = 0
    solutions = sp.solve(det_B, w2)
    print(f"\n\nSoluciones de la ecuación det(B) = 0: {solutions}")

    # Convertimos las soluciones simbólicas a números y calculamos la raíz cuadrada
    normal_frec = [sp.sqrt(sol.evalf()) for sol in solutions]
    normal_frec = np.array(normal_frec, dtype=np.complex_)

    print(f"\n\nLos modos normales son: {normal_frec} y {-1 * normal_frec}")

    return normal_frec

# Ejemplo de uso
Pendulos_super_simples = normal_frec(1, 1, 1, 1)