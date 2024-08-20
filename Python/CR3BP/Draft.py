import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Gravitational constant, m^3 kg^-1 s^-2
m1 = 5.972e24  # Mass of Earth, kg
m2 = 7.342e22  # Mass of Moon, kg
mu_earth = 398600.4418                   # km^3/s^2
mu_moon = 4902.800066                    # km^3/s^2
lu = 384400                              # Earth Moon distance km
mu = mu_moon / (mu_earth + mu_moon)      # Earth-Moon mass ratio
mu1 = 1 - mu
mu2 = mu

# Characteristic constants
tu = np.sqrt(lu**3 / (mu_earth + mu_moon))
vu = lu/tu
n = np.sqrt(G*(m1+m2)/lu**3)

# Equations of motion in the normalized CR3BP
def cr3bp_equations(t, state):
    x, y, z, dx, dy, dz = state
    t_s= [0, 10]
    t = np.linspace(t_s[0], t_s[1], 10000)
    r1 = np.sqrt((x + mu2*np.cos(n*t))**2 + (y + mu2*np.sin(n*t))**2 + z**2)
    r2 = np.sqrt((x - mu1*np.cos(n*t))**2 + (y - mu1*np.cos(n*t))**2 + z**2)
    
    ddx = - mu1 * (x + mu2*np.cos(n*t)) / r1**3 - mu2 * (x - mu1*np.cos(n*t)) / r2**3
    ddy = - mu1 * (y + mu2*np.sin(n*t)) / r1**3 - mu2 * (y - mu1*np.sin(n*t)) / r2**3
    ddz = - mu1 * z / r1**3 - mu2 * z / r2**3
    
    return [dx, dy, dz, ddx, ddy, ddz]

# Initial conditions (normalized)
ic = [
8.9833535487092597E-1,	-1.1587617145181475E-26,	1.7825269804266094E-35,	5.7536808647951106E-16,	4.7591168616820229E-1,	-1.2011348213495987E-34
]

# Time span for integration (Normalized)
t_span = [0, 10]
t_eval = np.linspace(t_span[0], t_span[1], 10000)

# Integrate the equations of motion
solution = solve_ivp(cr3bp_equations, t_span, ic, t_eval=t_eval, rtol=1e-10, atol=1e-12)

# Extract the trajectory
x = solution.y[0]*lu                     # Convert to km
y = solution.y[1]*lu                     # Convert to km
z = solution.y[2]*lu                     # Convert to km

# Plot the trajectory in 3D
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, label="Satellite Trajectory")
ax.plot([mu1*lu], [0], [0], 'ro', label="Moon")
ax.set_xlabel("x (km)")
ax.set_ylabel("y (km)")
ax.set_zlabel("z (km)")
ax.set_title("Moon Inertial Frame")
ax.legend()
plt.show()