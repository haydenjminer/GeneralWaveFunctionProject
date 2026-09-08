import math
import numpy as np
import matplotlib.pyplot as plt


# General Radial Wave Function Hydrogen
def laguerre(q,k, x):
    summation = 0
    for i in range(0, k+1):
        summation += ((-1) ** i) * ((math.factorial(k + q))/(math.factorial(k-i)*math.factorial(q+i))) * ((x**i)/ math.factorial(i))
    return summation


def general_radial_wave_function(x, a0, n, l, laguerre):
    term_1 = ((2/(n*a0))**3)
    term_2 = math.factorial(n-l-1)
    term_3 = 2*n*(math.factorial(n+l))
    term_4 = np.exp((-x)/(n*a0))
    term_5 =((2*x)/(n*a0))**l
    term_6 = laguerre(2*l+1, n-l-1, (2*x)/(n*a0))
    return math.sqrt(term_1*(term_2/term_3))*term_4*term_5*term_6


def radial_wave_plot(n, l, a0):
    x_variable = np.linspace(0, 70, 1000)
    x_output = x_variable ** 2
    function_output = general_radial_wave_function(x_variable, a0, n, l, laguerre)
    y_variable = 4 * math.pi * x_output * (abs(function_output)**2)
    plt.plot(x_variable, y_variable)
    plt.show()


# Spherical Harmonic Solutions 3d-Orbitals
# Spherical Harmonic 3d-z2 (Hydrogen-like)
def spherical_harmonic_3dz2():
    theta = np.linspace(0, 2 * np.pi, 100)
    r = ((3 * np.cos(theta) ** 2) - 1) * (np.sqrt(5 / (16 * np.pi)))
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    plt.polar(theta, abs(r))
    return plt.show()


# Spherical Harmonic 3d-xy (Hydrogen Like)
def spherical_harmonic_3dxy_plot_1():
    theta = np.linspace(0, 2*np.pi, 100)
    phi = np.pi/2
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    r = np.sqrt(15/(4*np.pi))*(np.sin(theta)**2)*np.cos(phi)*np.sin(phi)
    plt.polar(theta, r)
    return plt.show()


def spherical_harmonic_3dxy_plot_2():
    theta = np.linspace(0, 2*np.pi, 100)
    phi = 0
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    r = np.sqrt(15/(4*np.pi))*(np.sin(theta)**2)*np.cos(phi)*np.sin(phi)
    plt.polar(theta, r)
    return plt.show()


# Spherical Harmonic 3d-xz (Hydrogen Like)
def spherical_harmonic_3dxz_plot_1():
    theta = np.linspace(0, 2 * np.pi, 100)
    phi = (np.pi) / 2
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    r = np.sqrt(15/(4*np.pi)) * np.sin(theta) * np.cos(theta) * np.cos(phi)
    plt.polar(theta, abs(r))
    return plt.show()


def spherical_harmonic_3dxz_plot_2():
    theta = np.linspace(0, 2 * np.pi, 100)
    phi = 0
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    r = np.sqrt(15 / (4 * np.pi)) * np.sin(theta) * np.cos(theta) * np.cos(phi)
    plt.polar(theta, abs(r))
    return plt.show()


# Spherical Harmonic 3d-yz
def spherical_harmonic_3dyz_plot_1():
    theta = np.linspace(0, 2 * np.pi, 100)
    phi = (np.pi) / 2
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    r = np.sqrt(15/(4*np.pi)) * np.sin(theta) * np.cos(theta) * np.sin(phi)
    plt.polar(theta, abs(r))
    return plt.show()


def spherical_harmonic_3dyz_plot_2():
    theta = np.linspace(0, 2 * np.pi, 100)
    phi = 0
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    r = np.sqrt(15 / (4 * np.pi)) * np.sin(theta) * np.cos(theta) * np.sin(phi)
    plt.polar(theta, abs(r))
    return plt.show()


# Spherical Harmonic 3d_x2-y2
def spherical_harmonic_3dx2_y2_plot_1():
    theta = np.linspace(0, 2 * np.pi, 100)
    phi = np.pi/2
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    r = np.sqrt(15 / (16 * np.pi)) * (np.sin(theta)**2) * np.cos(2*phi)
    plt.polar(theta, r)
    return plt.show()

# Radial Wave Function Graphing
n_input = int(input("Enter quantum number n: "))
l_input = int(input("Enter quantum number l: "))
a0_defined = 0.529
print(radial_wave_plot(n_input, l_input, a0_defined))

graph_choice = input("Select a d-orbital graph to view: dz2, dxy, dxz, dyz, dx2-y2: ")
if graph_choice == "dz2":
    print(spherical_harmonic_3dz2())
elif graph_choice == "dxy":
    print(spherical_harmonic_3dxy_plot_1())
    print(spherical_harmonic_3dxy_plot_2())
elif graph_choice == "dxz":
    print(spherical_harmonic_3dxz_plot_1())
    print(spherical_harmonic_3dxz_plot_2())
elif graph_choice == "dyz":
    print(spherical_harmonic_3dyz_plot_1())
    print(spherical_harmonic_3dyz_plot_2())
elif graph_choice == "dx2-y2":
    print(spherical_harmonic_3dx2_y2_plot_1())
else:
    print("Invalid choice")