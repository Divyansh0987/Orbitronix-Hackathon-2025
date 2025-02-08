from Satellite import Satellite
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
from math import radians, cos, sin, asin, sqrt

# constraint = 1
# Best Parameters Found: {'target': 1024.0, 'params': {'arg_perigree': 218.55389831008577, 'inclination': -3.4565433999535458, 'raan': 293.8865442577971, 'semi_major_axis': 6745.052159307612}, 'constraint': 0.7492997252107478}
sat = Satellite(6745.052159307612, 0, -3.4565433999535458, 293.8865442577971, 218.55389831008577, 2 * np.pi * np.sqrt(6745.052159307612**3 / 3.986004418e5))

print("Orbital Period: ", 2 * np.pi * np.sqrt(6745.052159307612**3 / 3.986004418e5))

fig, ax = plt.subplots()

np.set_printoptions(threshold=sys.maxsize)

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r

min_distance = 999999
min_lat = 9999
min_long = 9999

for lat, long in zip(sat.get_ground_track(1, 10)[0], sat.get_ground_track(1, 10)[1]):
    test_distance = haversine(14, 80, lat, long)

    if test_distance < min_distance:
        min_distance = test_distance
        min_lat = lat
        min_long = long

print(min_distance)
print(min_lat)
print(min_long)

sat.plot_ground_track(sat.get_ground_track(1, 10), (1, 0, 0), fig, ax)
plt.show()