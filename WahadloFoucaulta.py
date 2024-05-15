from vpython import *

n = 15
g = vector(0, -9.8, 0)
L = 1
lat = pi / 2  # wychyl w bok
sol = 1000  # czas dnia

# Tworzenie sufitowego wspornika
ceiling = box(pos=vector(0, 1, 0), length=0.1, width=0.1, height=0.01, color=color.blue)
ball = sphere(pos=vector(0, 0, 0), radius=0.1, color=color.purple, make_trail=True, trail_type="points", interval=10, retain=50)
floor = cylinder(pos=vector(0, -1.1, 0), radius=3, axis=vector(0, -0.1, 0), color=color.green)
string = curve(radius=0.01)

markers = []
alpha = 0
while alpha < 2 * pi:
    marker = cylinder(pos=vector(0, -1, 0), axis=vector(cos(alpha), 0, sin(alpha)), radius=0.01, color=color.white)
    markers.append(marker)
    alpha += pi / 6

def pol2cart(theta):
    return vector(cos((3 * pi / 2 + theta)) * L, sin((3 * pi / 2 + theta)) * L, 0)

def latang(phi):
    return 2 * pi * sin(phi) / sol

# Tworzenie wykresu
graph_theta = graph(title="Pendulum Movement", xtitle="Time", ytitle="Position", fast=False)
massmovement = gcurve(graph=graph_theta, color=color.yellow)

L = 1
theta = pi / 9  # Ograniczenie wyhylu do 20 stopni
ball.pos = pol2cart(theta)
ball.m = 4  # masa
ball.L = vector(0, 0, 0)  # moment pedu
ball.I = ball.m * mag(ball.pos) ** 2  # moment bezwładnosci

string.append(pos=vector(0, 1, 0))
string.append(pos=ball.pos)

Fg = -ball.m * g
Tg = cross(ball.pos, Fg)

beta = 0
t = 0
dt = .01

scene.autoscale = False
scene.center = vector(0, 0, 0)
scene.range = 2
scene.camera.pos = vector(2, 1, 2)
scene.camera.axis = vector(-2, -1, -2)

light = distant_light(direction=vector(0.2, -0.2, 0.5), color=color.gray(0.8))

while True:
    rate(100)
    Tg = cross(ball.pos, Fg)
    ball.L = ball.L + Tg * dt
    theta = theta + ball.L.z / ball.I * dt

    ball.pos = pol2cart(theta)
    string.modify(1, ball.pos)

    if t < n:
        massmovement.plot(t, ball.pos.x * cos(beta))

    beta = beta + latang(lat)
    for marker in markers:
        marker.rotate(angle=latang(lat), axis=vector(0,1,0), origin=vector(0,1,0))

    t = t + dt
