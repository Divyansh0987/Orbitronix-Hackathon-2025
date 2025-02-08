import numpy as np
from bayes_opt import BayesianOptimization
from Satellite import Satellite
from scipy.optimize import NonlinearConstraint

R_Earth = 6371
mu_Earth = 3.986004418e5 # [km**3/s**2] Earth standard gravitational parameter

def get_optimal_close_approaches(semi_major_axis, inclination, raan, arg_perigree):
    # Define time in seconds
    totalTime = 2592000 # 1 month in seconds

    # Aqua orbit elements for close approach:
    aqua_inclination = 98.199
    aqua_raan = 95.206
    aqua_arg_perigree = 120.479
    aqua_semi_major_axis = R_Earth + 1

    # Compute orbital periods using Kepler's 3rd law
    T_orbit = 2 * np.pi * np.sqrt(semi_major_axis**3 / mu_Earth)
    aqua_T_orbit = 2 * np.pi * np.sqrt(aqua_semi_major_axis**3 / mu_Earth)
    n = 2 * np.pi / T_orbit
    aqua_n = 2 * np.pi / aqua_T_orbit

    i_rad = np.radians(inclination)
    raan_rad = np.radians(raan)
    w_rad = np.radians(arg_perigree)
    aqua_i_rad = np.radians(aqua_inclination)
    aqua_raan_rad = np.radians(aqua_raan)
    aqua_w_rad = np.radians(aqua_arg_perigree)

    # 2) TRUE ANOMALY & PERIFOCAL COORDINATES
    numPoints = totalTime
    t_vec = np.linspace(0, totalTime, numPoints)
    nu_vec = n * t_vec % (2 * np.pi)
    aqua_nu_vec = aqua_n * t_vec % (2 * np.pi)

    # Compute radius at each point
    r_pf = np.stack([semi_major_axis * np.cos(nu_vec), semi_major_axis * np.sin(nu_vec), np.zeros(numPoints)], axis=0)
    r_pf_aqua = np.stack([aqua_semi_major_axis * np.cos(aqua_nu_vec), aqua_semi_major_axis * np.sin(aqua_nu_vec), np.zeros(numPoints)], axis=0)

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

    idx = np.argwhere(np.diff(np.sign(np.degrees(theta_sph) - np.degrees(theta_sph_aqua)))).flatten()
    return len(t_vec[idx])

def constraint_function(semi_major_axis, inclination, raan, arg_perigree):
    T_orbit = 2 * np.pi * np.sqrt(semi_major_axis**3 / mu_Earth)
    sat = Satellite(semi_major_axis, 0, inclination, raan, arg_perigree, T_orbit)

    leastDist = 99999999

    for lat, long in zip(sat.get_ground_track(1, 10)[0], sat.get_ground_track(1, 10)[1]):
        _leastDist = np.sqrt(np.abs(14 - lat)**2 + np.abs(80 - long)**2)
        if _leastDist < leastDist:
            leastDist = _leastDist
    print(f"Latitude: {lat}, Longitude: {long}")

    return leastDist

constraint_limit = 1
constraint = NonlinearConstraint(constraint_function, -np.inf, constraint_limit)

pbounds = {
    'semi_major_axis': (R_Earth + 300, R_Earth + 700),  # Semi-major axis range
    'inclination': (-20, 20),  # Inclination range [0, 180]
    'raan': (0, 360),  # RAAN range [0, 360]
    'arg_perigree': (0, 360)  # Argument of Perigee range [0, 360]
}

# Run Bayesian Optimization
optimizer = BayesianOptimization(
    f=get_optimal_close_approaches,
    constraint=constraint,
    pbounds=pbounds,
    random_state=1,
)

optimizer.maximize(
    init_points=200,
    n_iter=20,
)

print("Best Parameters Found:", optimizer.max)