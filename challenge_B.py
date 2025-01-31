import numpy as np
from matplotlib import pyplot as plt

# 1) DEFINE ORBITAL PARAMETERS (ELLIPTICAL ORBIT FOR 1KM CLOSE APPROACH)
mu_Earth = 3.986004418e5  # [km**3/s**2] Earth standard gravitational parameter
R_Earth  = 6371           # [km] approximate Earth radius

# Define time in seconds
totalTime = 24 * 3600 # 24 hours in seconds

# Custom orbit elements for close approach:
incl_deg = 10        # Inclination i [deg]
raan_deg = 72        # Right Ascension of Ascending Node Ω [deg]
w_deg = 180          # Argument of Perigee ω [deg]
ecc = 0            # Eccentricity (elliptical orbit for close approach)
r_p = R_Earth + 1    # Perigee radius [km] (Earth radius + 1 km altitude)

# Aqua orbit elements for close approach:
aqua_incl_deg = 98.199
aqua_raan_deg = 95.206
aqua_w_deg = 120.479
aqua_ecc = 0
aqua_r_p = R_Earth + 1

# Compute semi-major axes
semi_major_axis = r_p / (1 - ecc)
aqua_semi_major_axis = aqua_r_p / (1 - aqua_ecc)

# Compute orbital periods using Kepler's 3rd law
T_orbit = 2 * np.pi * np.sqrt(semi_major_axis**3 / mu_Earth)
aqua_T_orbit = 2 * np.pi * np.sqrt(aqua_semi_major_axis**3 / mu_Earth)
n = 2 * np.pi / T_orbit
aqua_n = 2 * np.pi / aqua_T_orbit

def deg2rad(angle):
    return angle * np.pi / 180

# Convert angles to radians
i_rad = deg2rad(incl_deg)
raan_rad = deg2rad(raan_deg)
w_rad = deg2rad(w_deg)
aqua_i_rad = deg2rad(aqua_incl_deg)
aqua_raan_rad = deg2rad(aqua_raan_deg)
aqua_w_rad = deg2rad(aqua_w_deg)

# 2) TRUE ANOMALY & PERIFOCAL COORDINATES
numPoints = totalTime
t_vec = np.linspace(0, totalTime, numPoints)
nu_vec = n * t_vec % (2 * np.pi)
aqua_nu_vec = aqua_n * t_vec % (2 * np.pi)

# Compute radius at each point
r_orbit = semi_major_axis * (1 - ecc**2) / (1 + ecc * np.cos(nu_vec))
aqua_r_orbit = aqua_semi_major_axis * (1 - aqua_ecc**2) / (1 + aqua_ecc * np.cos(aqua_nu_vec))

# Perifocal frame position
r_pf = np.stack([r_orbit * np.cos(nu_vec), r_orbit * np.sin(nu_vec), np.zeros(numPoints)], axis=0)
r_pf_aqua = np.stack([aqua_r_orbit * np.cos(aqua_nu_vec), aqua_r_orbit * np.sin(aqua_nu_vec), np.zeros(numPoints)], axis=0)


# 3) ROTATE TO ECI FRAME
Rz_minusOmega = np.array([[np.cos(-raan_rad), np.sin(-raan_rad), 0], [-np.sin(-raan_rad), np.cos(-raan_rad), 0], [0, 0, 1]])

Rx_minusi = np.array([[1, 0, 0], [0, np.cos(-i_rad), np.sin(-i_rad)], [0, -np.sin(-i_rad), np.cos(-i_rad)]])

Rz_minusw = np.array([[np.cos(-w_rad), np.sin(-w_rad), 0], [-np.sin(-w_rad), np.cos(-w_rad), 0], [0, 0, 1]])

Q = Rz_minusOmega @ Rx_minusi @ Rz_minusw  # Use @ for matrix multiplication
r_eci = Q @ r_pf  # Use @ for matrix multiplication

Rz_minusOmega_aqua = np.array([[np.cos(-aqua_raan_rad), np.sin(-aqua_raan_rad), 0], [-np.sin(-aqua_raan_rad), np.cos(-aqua_raan_rad), 0], [0, 0, 1]])

Rx_minusi_aqua = np.array([[1, 0, 0], [0, np.cos(-aqua_i_rad), np.sin(-aqua_i_rad)], [0, -np.sin(-aqua_i_rad), np.cos(-aqua_i_rad)]])

Rz_minusw_aqua = np.array([[np.cos(-aqua_w_rad), np.sin(-aqua_w_rad), 0], [-np.sin(-aqua_w_rad), np.cos(-aqua_w_rad), 0], [0, 0, 1]])

Qa = Rz_minusOmega_aqua @ Rx_minusi_aqua @ Rz_minusw_aqua  # Use @ for matrix multiplication
r_eci_aqua = Qa @ r_pf_aqua  # Use @ for matrix multiplication

# Extract coordinates (no change needed here)
x_eci = r_eci[0, :]
y_eci = r_eci[1, :]
z_eci = r_eci[2, :]

x_eci_aqua = r_eci_aqua[0, :]
y_eci_aqua = r_eci_aqua[1, :]
z_eci_aqua = r_eci_aqua[2, :]

# 4) CONVERT ECI TO SPHERICAL COORDINATES (no change needed here)
r_sph = np.sqrt(x_eci**2 + y_eci**2 + z_eci**2)
theta_sph = np.arctan2(np.sqrt(x_eci**2 + y_eci**2), z_eci)
phi_sph = np.arctan2(y_eci, x_eci)

r_sph_aqua = np.sqrt(x_eci_aqua**2 + y_eci_aqua**2 + z_eci_aqua**2)
theta_sph_aqua = np.arctan2(np.sqrt(x_eci_aqua**2 + y_eci_aqua**2), z_eci_aqua)
phi_sph_aqua = np.arctan2(y_eci_aqua, x_eci_aqua)

# 5) PLOTTING
fig, ax = plt.subplots(2, 1, figsize=(8, 6))
fig.suptitle('Spherical Coordinates', fontsize=12, color='w')

# Plot the polar angle
ax[0].plot(t_vec, np.degrees(theta_sph), 'r', linewidth=1.5, label='Data 1') #Added label
ax[0].plot(t_vec, np.degrees(theta_sph_aqua), 'b', linewidth=1.5, label='Data 2') #Added label
ax[0].set_xlabel('Time [s]')
ax[0].set_ylabel(r'$theta$ [deg]') #Use raw string for LaTeX formatting
ax[0].set_title('Polar Angle vs. Time')
ax[0].grid(True)
ax[0].legend() #Show legend for the plot

# Plot the azimuth angle
ax[1].plot(t_vec, np.degrees(phi_sph), 'r', linewidth=1.5, label='Data 1')  #Added label
ax[1].plot(t_vec, np.degrees(phi_sph_aqua), 'b', linewidth=1.5, label='Data 2')  #Added label
ax[1].set_xlabel('Time [s]')
ax[1].set_ylabel(r'$phi$ [deg]') #Use raw string for LaTeX formatting
ax[1].set_title('Azimuth Angle vs. Time')
ax[1].grid(True)
ax[1].legend() #Show legend for the plot

idx = np.argwhere(np.diff(np.sign(np.degrees(theta_sph) - np.degrees(theta_sph_aqua)))).flatten()
# plt.plot(t_vec[idx], np.sign(np.degrees(theta_sph))[idx], 'ro')

print(len(t_vec[idx]))

# Adjust spacing between subplots
plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # adjust to prevent overlap with suptitle

# Show the plot
plt.show()