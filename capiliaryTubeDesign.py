import math
import cmath
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# Input

Tc = 40  # Condensing temperature in
Te = 18  # Evaporator temperature
m = 0.0256  # Mass Flow rate Kg/s
d = 0.094  # Diameter in inches
D = 0.0254*d
e = (1e-6)*0.032
A = 3.14*(D*D)/4  # Area in Mtere Square
G = m/A

# Created array for each variable quantity to get graph and discretized value
reynolds = []
specificVolume = []
specificEnthalpy = []
viscosity = []
pressure = []
velocity = []
o = []
n = []


def EqRoots(a, b, c):

    dis = b*b-4*a*c
    sqrt_val = np.sqrt(dis)

    # checking condition for discriminant
    if dis > 0:

        r1 = (-b + sqrt_val)/(2 * a)
        r2 = (-b - sqrt_val)/(2 * a)
        if r1 > 0:
            return r1
        if r2 > 0:
            return r2
        else:
            return 0

    elif dis == 0:

        return (-b / (2 * a))


df = pd.read_csv("R410A.csv")  # File for R410A location should be in same folder 
df.set_index('Temp')

prop = df.loc[df["Temp"] == Tc]
arr = prop.values.tolist()

p1 = 1e6*arr[0][1]
h1 = 1e3*arr[0][5]
v1 = arr[0][3]
u1 = arr[0][11]
Cpl = 1e3*arr[0][9]

V1 = G*v1
Re1 = (V1*D/(u1*v1))
f1 = 0.25*(math.log((math.e / (3.7*D))+(5.74/(Re1 ** 0.9)))) ** (-2)

for i in range(Tc, Te-1, -1):
    o.append(i)


for j in range(len(o)):
    prop2 = df.loc[df["Temp"] == o[j]]
    arr2 = prop2.values.tolist()
    p2 = 1e6*arr2[0][1]
    hf2 = 1e3*arr2[0][5]
    hg2 = 1e3*arr2[0][6]
    vf2 = arr2[0][3]
    vg2 = arr2[0][4]
    uf2 = arr2[0][11]
    ug2 = arr2[0][12]

    a = 0.5*(vg2-vf2) * (vg2-vf2)*G * G
    b = (hg2-hf2)+(vf2*(vg2-vf2))*G * G
    c = (hf2-h1)+0.5*(vf2 * vf2*G * G)-(V1 * V1/2)

    x2 = EqRoots(a, b, c)

    # print(f"x2 is {x2}")

    v2 = vf2*(1-x2) + vg2*x2
    h2 = hf2*(1-x2) + hg2*x2
    u2 = uf2*(1-x2) + ug2*x2
    V2 = G*v2

    Re2 = (V2*D/(u2*v2))
    f2 = 0.25*(math.log((math.e / (3.7*D))+(5.74/(Re2 ** 0.9)))) ** (-2)

    Vm = (V1+V2)*0.5
    fm = (f1+f2)*0.5

    # print("Vm", Vm, "\n", "fm", fm,)

    L12 = ((2)*D*((p1-p2)-G*(V2-V1))/(fm*Vm*G))
    n.append(L12)
    velocity.append(Vm)
    pressure.append(p2)
    reynolds.append(Re2)
    specificVolume.append(v2)
    specificEnthalpy.append(h2)
    viscosity.append(u2)


# print(o[j])

print(" Lenght of Capillary Tube in meter is ", L12)

# print(n)

plt.plot(n, o)
plt.xlabel('Capillary tube length (m)')
plt.ylabel('temperature (deg C)')
plt.title('temperature vs Capillary tube length')
plt.show()

plt.plot(n, velocity)
plt.xlabel('Capillary tube length (m)')
plt.ylabel('velocity (m/s)')
plt.title('Velocity vs Capillary tube length')
plt.show()

plt.plot(n, pressure)
plt.xlabel('Capillary tube length (m)')
plt.ylabel('pressure')
plt.title(' pressure vs Capillary tube length')
plt.show()

plt.plot(n, reynolds)
plt.xlabel('Capillary tube length (m)')
plt.ylabel('Reynolds number')
plt.title('Reynolds number vs Capillary tube length')
plt.show()

plt.plot(n, specificVolume)
plt.xlabel('Capillary tube length (m)')
plt.ylabel('specific volume')
plt.title('specific volume vs Capillary tube length')
plt.show()

plt.plot(n, specificEnthalpy)
plt.xlabel('Capillary tube length (m)')
plt.ylabel(' specific enthalpy')
plt.title('specific enthalpy vs Capillary tube length')
plt.show()

plt.plot(n, viscosity)
plt.xlabel('Capillary tube length (m)')
plt.ylabel('dynamic viscosity in (Pa.s)')
plt.title('Dynamic viscosity vs Capillary tube length')
plt.show()
