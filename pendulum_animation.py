#zeby odpalic w terminalu : python pendulum_animation.py

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.animation import FuncAnimation

# Parametry fizyczne
L = 0.3  # Długość wahadła (m)
g = 9.81  # Przyspieszenie grawitacyjne (m/s^2)
T = 10  # Okres obrotu Ziemi (s)
Ω = 2 * np.pi / T  # Pulsacja obrotu Ziemi
l = np.radians(45)  # Szerokość geograficzna wahadła
ω0 = np.sqrt(g / L)  # Pulsacja własna wahadła

# Funkcja rozwiązująca równania różniczkowe za pomocą odeint
def solve_ODE(equation, f0, t_values, colonne):
    solution = odeint(equation, f0, t_values)
    f_values = solution[:, colonne]
    return f_values

# Funkcja reprezentująca równania różniczkowe wahadła.
def Eq_DP(y, t):
    θ, dθdt, φ, dφdt = y

    # Równania różniczkowe dla θ i φ
    d2θdt2 = 2 * Ω * np.sin(l) * np.sin(θ) * np.cos(θ) * dφdt - 2 * Ω * np.sin(φ) * np.sin(θ) ** 2 * np.cos(
        l) * dφdt - ω0 ** 2 * np.sin(θ) + np.sin(θ) * np.cos(θ) * dφdt ** 2

    d2φdt2 = (-2 * Ω * np.sin(l) * np.cos(θ) * dθdt + 2 * Ω * np.sin(φ) * np.sin(θ) * np.cos(
        l) * dθdt - 2 * np.cos(θ) * dθdt * dφdt) / np.sin(θ)

    return [dθdt, d2θdt2, dφdt, d2φdt2]

# Warunki początkowe
θ_0 = np.radians(60)
φ_0 = np.radians(0.0)
dθdt_0 = np.radians(0.0)
dφdt_0 = np.radians(0.0)
initial_conditions = [θ_0, φ_0, dθdt_0, dφdt_0]

# Przedział czasowy dla symulacji
t_start = 0.0
t_end = 60
t_values = np.arange(t_start, t_end, 0.02)

# Rozwiązanie równań różniczkowych dla wahadła.
θ = solve_ODE(Eq_DP, initial_conditions, t_values, 0)
φ = solve_ODE(Eq_DP, initial_conditions, t_values, 2)

# Współrzędne kartezjańskie wahadła
x = L * np.sin(θ) * np.cos(φ)
y = L * np.sin(θ) * np.sin(φ)
z = -L * np.cos(θ)

# Liczba mas dodatkowych na linie wahadła
n = 100

# Utworzenie figury 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([-L, L])
ax.set_ylim([-L, L])
ax.set_zlim([-L, 0])

# Inicjalizacja linii trajektorii, kulki wahadła i liny
trajectory_line, = ax.plot([], [], [], color='blue')
pendulum_bob, = ax.plot([], [], [], color='red', marker='o', markersize=10)
rope = [ax.plot([], [], [], color='black', marker='o', markersize=1)[0] for _ in range(n)]

# Funkcja animacji
def animate(i):
    trajectory_line.set_data(x[:i], y[:i])
    trajectory_line.set_3d_properties(z[:i])
    pendulum_bob.set_data([x[i]], [y[i]])
    pendulum_bob.set_3d_properties([z[i]])

    # Aktualizacja pozycji dodatkowych mas na linie
    for j in range(n):
        # Obliczenie pozycji
        x_additional = (L - j * L / n) * np.sin(θ[i]) * np.cos(φ[i])
        y_additional = (L - j * L / n) * np.sin(θ[i]) * np.sin(φ[i])
        z_additional = -(L - j * L / n) * np.cos(θ[i])

        # Aktualizacja pozycji dodatkowej masy
        rope[j].set_data([x_additional], [y_additional])
        rope[j].set_3d_properties([z_additional])

    return [trajectory_line, pendulum_bob] + rope

# Animacja trajektorii
ani = FuncAnimation(fig, animate, frames=len(t_values), interval=1, blit=True)

plt.show()
