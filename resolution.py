import math
import numpy as np
import matplotlib.pyplot as plt
import pdb


plt.ion()

def RayleighCriteria(d: float, wavelength: float=500e-9) -> float:
    return (1.22*wavelength/d)*(180.0/math.pi)*3600.0


def seeingDiskToDiameter(d: float) -> float:
    return 1.1902439024390246*d

diameters = np.logspace(-3, 1, 1000)
resolutions = np.array([RayleighCriteria(d) for d in diameters])

eye_diameters = np.array([5.0e-3, 9.0e-3])
eye_resolutions = np.array([5.0*1e-4*(180.0/math.pi)*3600.0, 2.0*1e-4*(180.0/math.pi)*3600.0])

plt.figure(0)
plt.xlabel('diameter (mm)')
plt.ylabel('resolution (arcsec)')
plt.loglog(1000*diameters, resolutions, label='Rayleigh criteria resolution')
plt.scatter(1000*eye_diameters, eye_resolutions, label='Human eye')
plt.hlines()
plt.show()

pdb.set_trace()