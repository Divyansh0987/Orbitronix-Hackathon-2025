from Satellite import Satellite
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('agg')

PDFs = []

sat_g = Satellite(7080.6, 0, 98.199, 95.206, 0, 5940)
sat_h = Satellite(7080.6, 0, 34, 95.206, 0, 5940)

print("Making GSAT plots...")

fig, ax = plt.subplots()
sat_g.plot_ground_track(sat_g.get_ground_track(1, 10), (1, 0, 0), fig, ax)
plt.savefig("GSAT-ground-track.pdf")

fig, ax = plt.subplots()
sat_g.plot_coverage_area(sat_g.get_ground_track(1, 10), (1, 0, 0), fig, ax)
plt.savefig("GSAT-coverage-area.pdf")

print("GSAT plots done.")

print("Making HSAT plots...")

fig, ax = plt.subplots()
sat_h.plot_ground_track(sat_h.get_ground_track(1, 10), (1, 0, 0), fig, ax)
plt.savefig("HSAT-ground-track.pdf")

fig, ax = plt.subplots()
sat_h.plot_coverage_area(sat_h.get_ground_track(1, 10), (1, 0, 0), fig, ax)
plt.savefig("HSAT-coverage-area.pdf")

print("HSAT plots done.")

print("All done")

plt.close()
