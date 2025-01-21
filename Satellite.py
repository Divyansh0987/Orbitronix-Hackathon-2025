import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Circle
from matplotlib.colors import colorConverter

class Satellite:
    earth_radius = 6371
    earth_rot_speed = 360 / 86164

    def __init__(self, semi_major_axis, e, inclination, raan,  arg_perigree, period):
        self.semi_major_axis = semi_major_axis
        self.e = e
        self.inclination = inclination
        self.raan = raan
        self.arg_perigree = arg_perigree
        self.period = period
        self.mean_motion = 360 / self.period
    
    def get_ground_track(self, num_orbits, time_step):
        latitudes = []
        longitudes = []
        self.num_orbits = num_orbits
        self.time_step = time_step
        time_array = np.arange(0, num_orbits * self.period, time_step)
        for time in time_array:
            mean_anomaly = self.mean_motion * time

            nu = np.radians(mean_anomaly % 360)

            x_orb = self.semi_major_axis * np.cos(nu)
            y_orb = self.semi_major_axis * np.sin(nu)
            z_orb = 0

            X_eci = (x_orb * (np.cos(self.raan) * np.cos(Satellite.earth_rot_speed * time) - np.sin(self.raan) * np.sin(Satellite.earth_rot_speed * time) * np.cos(self.inclination))
                 - y_orb * np.sin(Satellite.earth_rot_speed * time))
            Y_eci = (x_orb * (np.sin(self.raan) * np.cos(Satellite.earth_rot_speed * time) + np.cos(self.raan) * np.sin(Satellite.earth_rot_speed * time) * np.cos(self.inclination))
                    + y_orb * np.cos(Satellite.earth_rot_speed * time))
            Z_eci = x_orb * np.sin(self.inclination)

            r = np.sqrt(X_eci**2 + Y_eci**2 + Z_eci**2)
            latitude = np.degrees(np.arcsin(Z_eci / r))
            longitude = (np.degrees(np.arctan2(Y_eci, X_eci)) - Satellite.earth_rot_speed * time) % 360
            
            if longitude > 180:
                longitude -= 360
            
            latitudes.append(latitude)
            longitudes.append(longitude)
        return np.array(latitudes), np.array(longitudes)

    def get_coverage_radius(self):
        altitude = self.semi_major_axis - Satellite.earth_radius
        alpha = np.arcsin(Satellite.earth_radius / (Satellite.earth_radius + altitude))
        radius = alpha * altitude / (90 - alpha)

        return radius

    def plot_ground_track(self, ground_track, color, fig, ax):
        latitudes, longitudes = ground_track

        m = Basemap(projection='cyl', resolution='l', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)
        m.drawcoastlines()
        m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180, 181, 60), labels=[0, 0, 0, 1])

        # Plot and connect ground track points
        x, y = m(longitudes, latitudes)
        ax.plot(x, y, marker='o', markersize=1, linewidth=0.5, color=color, label='Ground Track')

    def plot_coverage_area(self, ground_track, color, fig, ax):
        latitudes, longitudes = ground_track

        m = Basemap(projection='cyl', resolution='l', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)
        m.drawcoastlines()
        m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180, 181, 60), labels=[0, 0, 0, 1])

        # Highlight the coverage area for one orbit
        for latitude in latitudes:
            for longitude in longitudes:
                circle = plt.Circle(m(latitude, longitude), self.get_coverage_radius(), color=color, fill=True, alpha = 0.1)
                ax.add_artist(circle)