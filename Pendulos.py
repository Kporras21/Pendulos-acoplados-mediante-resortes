import numpy as np
import sympy as sp
from scipy.constants import g
import cmath
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

def normal(k1, k2, m, L):
    """
        Función que recibe la masa colgante de ambos pendulos, la longitud L de la cuerda de los pendulos y las constantes de contracción
    de los resortes que acoplan el sistema y devuelve; las frecuencias características, los modos normales y su correlación 
    además de un grafico con el movimiento de ambos pendulos.

    Args:
        k1 (float): Constante del resorte pendulo pared.
        k2 (float): Constante del resorte pendulo pendulo.
        m (float): Masa de los pendulos.
        L (float): Longitud de la cuerda de cada pendulo.

    Return:
        times (array) = arreglo de tiempos donde se calculó númericamente el movimiento.
        theta_11 (array): Arreglo de posiciones para el péndulo de la izquierda en su primer modo.
        theta_12 (array): Arreglo de posiciones para el péndulo de la derecha en su primer modo.
        theta_21 (array): Arreglo de posiciones para el péndulo de la izquierda en su segundo modo modo.
        theta_22 (array): Arreglo de posiciones para el péndulo de la derecha en su segundo modo.








    """
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


    # Definimos las equaciones con las frecuencias normales al cuadrado
    a11 = (m * g * L + (L**2) * k2 - solutions[0] * m * (L**2))/(-k2 * (L**2))


    print(f'La relación entre los modos para la primera frecuencia es de: a1 = {a11} *a2\n')

    a12 = (m * g * L + (L**2) * k2 - solutions[1] * m * (L**2))/(-k2 * (L**2))

    print(f'La relación entre los modos para la primera frecuencia es de: a1 = {a12} * a2\n\n')
    


    return normal_frec, a11, a12

frecuencias_U , A11, A12 = normal(1,1,1,1)

def movement(freqs, mode_1, mode_2, A, t_max):
    times = np.linspace(0,t_max, 20)
    f1 = freqs[0].real
    f2 = freqs[1].real

    theta_11 = np.zeros(np.size(times))
    theta_21 = np.zeros(np.size(times))
    
    for i in range(len(times)-1):
        theta_11[i] = A * mode_1 * (np.cos((f1)* times[i]))
        theta_21[i] = A * (np.cos((f1)*times[i]))
    

    plt.plot(times, theta_11, color = 'blue')
    plt.plot(times, theta_21, color = 'red')
    plt.grid()
    plt.show()
    

    theta_12 = np.zeros(np.size(times))
    theta_22 = np.zeros(np.size(times))

    for i in range(len(times)-1):
        theta_12[i] = A * mode_2 * (np.cos((f2)*times[i]))
        theta_22[i] = A * (np.cos((f2)*times[i]))
    plt.plot(times, theta_12, color = 'blue')
    plt.plot(times, theta_22, color = 'red')
    plt.grid()
    plt.show()

    return times, theta_11, theta_21, theta_12, theta_22

times, theta_11, theta_21, theta_12, theta_22 = movement(frecuencias_U, A11, A12, 5, 30)




