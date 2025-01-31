from Satellite import Satellite
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# constraint = 1
# Best Parameters Found: {'target': 1024.0, 'params': {'arg_perigree': 218.55389831008577, 'inclination': -3.4565433999535458, 'raan': 293.8865442577971, 'semi_major_axis': 6745.052159307612}, 'constraint': 0.7492997252107478}
sat = Satellite(6745.052159307612, 0, -3.4565433999535458, 293.8865442577971, 218.55389831008577, 2 * np.pi * np.sqrt(6745.052159307612**3 / 3.986004418e5))

print("Orbital Period: ", 2 * np.pi * np.sqrt(6745.052159307612**3 / 3.986004418e5))

fig, ax = plt.subplots()
sat.plot_ground_track(sat.get_ground_track(1, 10), (1, 0, 0), fig, ax)
plt.show()