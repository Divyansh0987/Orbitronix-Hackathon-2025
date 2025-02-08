import numpy as np
from bayes_opt import BayesianOptimization

seed = 42
np.random.seed(seed)

# Define global Earth radius (fixes NameError)
R_Earth = 6371  # km, approximate Earth radius

def get_optimal_close_approaches(semi_major_axis, inclination, raan, arg_perigee):
    mu_Earth = 3.986004418e5  # [km^3/s^2] Earth standard gravitational parameter

    # Define total time for simulation (e.g., 30 days)
    totalTime = 30 * 24 * 3600  # 30 days in seconds

    # Aqua satellite orbit parameters
    aqua_inclination = 98.199
    aqua_raan = 95.206
    aqua_arg_perigee = 120.479
    aqua_semi_major_axis = R_Earth + 1  # Fix: Uses global R_Earth

    # Compute orbital periods using Kepler’s 3rd law
    T_orbit = 2 * np.pi * np.sqrt(semi_major_axis**3 / mu_Earth)
    aqua_T_orbit = 2 * np.pi * np.sqrt(aqua_semi_major_axis**3 / mu_Earth)
    n = 2 * np.pi / T_orbit
    aqua_n = 2 * np.pi / aqua_T_orbit

    # Convert angles to radians
    i_rad = np.radians(inclination)
    raan_rad = np.radians(raan)
    w_rad = np.radians(arg_perigee)
    aqua_i_rad = np.radians(aqua_inclination)
    aqua_raan_rad = np.radians(aqua_raan)
    aqua_w_rad = np.radians(aqua_arg_perigee)

    # Define time vector
    numPoints = 10000  # Reduce computation by using fewer points
    t_vec = np.linspace(0, totalTime, numPoints)
    nu_vec = (n * t_vec) % (2 * np.pi)
    aqua_nu_vec = (aqua_n * t_vec) % (2 * np.pi)

    # Compute orbital positions (assuming circular orbit)
    r_pf = np.stack([semi_major_axis * np.cos(nu_vec), semi_major_axis * np.sin(nu_vec), np.zeros(numPoints)], axis=0)
    r_pf_aqua = np.stack([aqua_semi_major_axis * np.cos(aqua_nu_vec), aqua_semi_major_axis * np.sin(aqua_nu_vec), np.zeros(numPoints)], axis=0)

    # Compute rotation matrices
    Rz_minusOmega = np.array([[np.cos(-raan_rad), np.sin(-raan_rad), 0], [-np.sin(-raan_rad), np.cos(-raan_rad), 0], [0, 0, 1]])
    Rx_minusi = np.array([[1, 0, 0], [0, np.cos(-i_rad), np.sin(-i_rad)], [0, -np.sin(-i_rad), np.cos(-i_rad)]])
    Rz_minusw = np.array([[np.cos(-w_rad), np.sin(-w_rad), 0], [-np.sin(-w_rad), np.cos(-w_rad), 0], [0, 0, 1]])

    Q = Rz_minusOmega @ Rx_minusi @ Rz_minusw
    r_eci = Q @ r_pf

    Rz_minusOmega_aqua = np.array([[np.cos(-aqua_raan_rad), np.sin(-aqua_raan_rad), 0], [-np.sin(-aqua_raan_rad), np.cos(-aqua_raan_rad), 0], [0, 0, 1]])
    Rx_minusi_aqua = np.array([[1, 0, 0], [0, np.cos(-aqua_i_rad), np.sin(-aqua_i_rad)], [0, -np.sin(-aqua_i_rad), np.cos(-aqua_i_rad)]])
    Rz_minusw_aqua = np.array([[np.cos(-aqua_w_rad), np.sin(-aqua_w_rad), 0], [-np.sin(-aqua_w_rad), np.cos(-aqua_w_rad), 0], [0, 0, 1]])

    Qa = Rz_minusOmega_aqua @ Rx_minusi_aqua @ Rz_minusw_aqua
    r_eci_aqua = Qa @ r_pf_aqua

    # Compute spherical coordinates
    theta_sph = np.arctan2(np.sqrt(r_eci[0]**2 + r_eci[1]**2), r_eci[2])
    theta_sph_aqua = np.arctan2(np.sqrt(r_eci_aqua[0]**2 + r_eci_aqua[1]**2), r_eci_aqua[2])

    # Count close approaches
    close_approaches = np.sum(np.abs(np.degrees(theta_sph) - np.degrees(theta_sph_aqua)) < 1)

    return close_approaches

# Define the bounds for Bayesian Optimization (R_Earth is now defined globally)
pbounds = {
    'semi_major_axis': (R_Earth + 1, R_Earth + 500),  # Semi-major axis range
    'inclination': (0, 180),  # Inclination range [0, 180]
    'raan': (0, 360),  # RAAN range [0, 360]
    'arg_perigee': (0, 360)  # Argument of Perigee range [0, 360]
}

# Run Bayesian Optimization
optimizer = BayesianOptimization(
    f=get_optimal_close_approaches,
    pbounds=pbounds,
    random_state=1,
)

optimizer.maximize(
    init_points=2000,
    n_iter=20,
)

print("Best Parameters Found:", optimizer.max)
