# Capillary Tube Sizing Script for R410 Refrigerant
# Description: Calculates required capillary tube length for a given mass flow rate and operating conditions.

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# -------------------------- Input Parameters --------------------------
Tc = 40  # Condensing temperature (°C)
Te = 18  # Evaporator temperature (°C)
m = 0.0256  # Mass flow rate (kg/s)
d_inch = 0.094  # Diameter in inches
D = 0.0254 * d_inch  # Diameter in meters
A = math.pi * (D**2) / 4  # Cross-sectional area in m^2
G = m / A  # Mass flux

# -------------------------- Data Initialization --------------------------
tube_lengths, velocity, pressure = [], [], []
reynolds, specific_volume, specific_enthalpy = [], [], []
viscosity, temps = [], []

# -------------------------- Helper Function --------------------------
def solve_quadratic_root(a, b, c):
    dis = b**2 - 4*a*c
    if dis < 0:
        return 0
    sqrt_val = np.sqrt(dis)
    r1 = (-b + sqrt_val) / (2*a)
    r2 = (-b - sqrt_val) / (2*a)
    if r1 > 0:
        return r1
    elif r2 > 0:
        return r2
    else:
        return 0

# -------------------------- Load R410 Data --------------------------
csv_filename = "R410.csv"

if not os.path.exists(csv_filename):
    print(f"'{csv_filename}' not found. Creating sample dataset...")
    sample_data = {
        "Temp": list(range(10, 51)),
        "p_vap (MPa)": np.linspace(0.5, 2.6, 41),
        "hf": np.linspace(100, 300, 41),
        "hg (kJ/kg)": np.linspace(400, 600, 41),
        "vf (m3/kg)": np.linspace(0.0008, 0.0012, 41),
        "vg (m3/kg)": np.linspace(0.02, 0.04, 41),
        "muf (Pa-s)": np.linspace(1.2e-5, 1.8e-5, 41),
        "mug (Pa-s)": np.linspace(9e-6, 1.1e-5, 41)
    }
    df = pd.DataFrame(sample_data)
    df.to_csv(csv_filename, index=False)
else:
    df = pd.read_csv(csv_filename)

df["Temp"] = df["Temp"].astype(int)

try:
    row_init = df[df["Temp"] == Tc].iloc[0]
except IndexError:
    raise ValueError(f"Temperature {Tc}°C not found in {csv_filename}. Check your data.")

p1 = row_init["p_vap (MPa)"] * 1e6
h1 = row_init["hf"] * 1e3
v1 = row_init["vf (m3/kg)"]
u1 = row_init["muf (Pa-s)"]

V1 = G * v1
Re1 = (V1 * D) / (u1 * v1)
f1 = 0.25 / (math.log((math.e / (3.7 * D)) + (5.74 / Re1**0.9)))**2

# -------------------------- Main Calculation Loop --------------------------
length = 0
for T in range(Tc, Te - 1, -1):
    row = df[df["Temp"] == T]
    if row.empty:
        continue
    row = row.iloc[0]

    p2 = row["p_vap (MPa)"] * 1e6
    hf2 = row["hf"] * 1e3
    hg2 = row["hg (kJ/kg)"] * 1e3
    vf2 = row["vf (m3/kg)"]
    vg2 = row["vg (m3/kg)"]
    uf2 = row["muf (Pa-s)"]
    ug2 = row["mug (Pa-s)"]

    a = 0.5 * (vg2 - vf2)**2 * G**2
    b = (hg2 - hf2) + vf2 * (vg2 - vf2) * G**2
    c = (hf2 - h1) + 0.5 * (vf2**2 * G**2) - (V1**2 / 2)

    x2 = solve_quadratic_root(a, b, c)
    v2 = vf2 * (1 - x2) + vg2 * x2
    h2 = hf2 * (1 - x2) + hg2 * x2
    u2 = uf2 * (1 - x2) + ug2 * x2
    V2 = G * v2

    Re2 = (V2 * D) / (u2 * v2)
    f2 = 0.25 / (math.log((math.e / (3.7 * D)) + (5.74 / Re2**0.9)))**2

    Vm = 0.5 * (V1 + V2)
    fm = 0.5 * (f1 + f2)

    dL = (2 * D * ((p1 - p2) - G * (V2 - V1))) / (fm * Vm * G)
    length += dL

    tube_lengths.append(length)
    velocity.append(Vm)
    pressure.append(p2)
    reynolds.append(Re2)
    specific_volume.append(v2)
    specific_enthalpy.append(h2)
    viscosity.append(u2)
    temps.append(T)

# -------------------------- Final Output --------------------------
if tube_lengths:
    print(f"\nEstimated Capillary Tube Length: {tube_lengths[-1]:.4f} m\n")
else:
    print("\nNo valid temperature data found. Ensure R410.csv contains required temperatures.\n")

# -------------------------- Plotting --------------------------
def plot_graph(x, y, xlabel, ylabel, title):
    plt.figure()
    plt.plot(x, y, marker='o')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if tube_lengths:
    plot_graph(temps, tube_lengths, 'Temperature (°C)', 'Capillary tube length (m)', 'Capillary Tube Length vs Temperature')
    plot_graph(tube_lengths, velocity, 'Capillary tube length (m)', 'Velocity (m/s)', 'Velocity vs Capillary Tube Length')
    plot_graph(tube_lengths, pressure, 'Capillary tube length (m)', 'Pressure (Pa)', 'Pressure vs Capillary Tube Length')
    plot_graph(tube_lengths, reynolds, 'Capillary tube length (m)', 'Reynolds Number', 'Reynolds Number vs Capillary Tube Length')
    plot_graph(tube_lengths, specific_volume, 'Capillary tube length (m)', 'Specific Volume (m3/kg)', 'Specific Volume vs Capillary Tube Length')
    plot_graph(tube_lengths, specific_enthalpy, 'Capillary tube length (m)', 'Specific Enthalpy (J/kg)', 'Specific Enthalpy vs Capillary Tube Length')
    plot_graph(tube_lengths, viscosity, 'Capillary tube length (m)', 'Dynamic Viscosity (Pa.s)', 'Dynamic Viscosity vs Capillary Tube Length')
