import numpy as np
import sympy as sp
from scipy.constants import g
import cmath
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

def normal(k1, k2, m, L):
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
    print(theta_11)

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

# Crear la figura y los ejes
fig, ax = plt.subplots(2, 1, figsize=(8, 8))

# Inicializar las líneas que se actualizarán en la animación
line1, = ax[0].plot([], [], lw=2, color='blue', label='Modo 1')
line2, = ax[0].plot([], [], lw=2, color='red', label='Modo 2')
line3, = ax[1].plot([], [], lw=2, color='blue', label='Modo 1')
line4, = ax[1].plot([], [], lw=2, color='red', label='Modo 2')

# Configurar los ejes
ax[0].set_xlim(0, 30)
ax[0].set_ylim(-5, 5)
ax[0].set_title('Primera Frecuencia')
ax[0].legend()
ax[0].grid()

ax[1].set_xlim(0, 30)
ax[1].set_ylim(-5, 5)
ax[1].set_title('Segunda Frecuencia')
ax[1].legend()
ax[1].grid()

# Función de inicialización
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    return line1, line2, line3, line4

# Función de actualización
def update(frame):
    line1.set_data(times[:frame], theta_11[:frame])
    line2.set_data(times[:frame], theta_21[:frame])
    line3.set_data(times[:frame], theta_12[:frame])
    line4.set_data(times[:frame], theta_22[:frame])
    return line1, line2, line3, line4

# Crear la animación
ani = animation.FuncAnimation(fig, update, frames=len(times), init_func=init, blit=True)

# Guardar la animación como un archivo mp4
#ani.save('pendulum_animation.mp4', writer='ffmpeg', fps=30)

# Mostrar la animación
plt.tight_layout()
plt.show()

