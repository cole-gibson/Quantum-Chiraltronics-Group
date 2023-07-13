import numpy as np
import os
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300

from scipy.spatial import Delaunay
sublattices = 1 # set to 1 if no sublattices

atominfo = np.loadtxt("atoms-coords.data", usecols=range(5), skiprows=1, dtype=np.float32)
allpoints = atominfo[:,2:5]
# Magnetic information about individual spins - takes vectors in form (S_x,S_y,S_z)
allspins = np.loadtxt("spins-00000030.data", usecols=range(3), skiprows=1, dtype=np.float32)
# The above files need to be arranged so that row n of the spin array (ie the nth spin vector) is the
# spin corresponding to the lattice point given by row n of the points array. This is what allows us to use the
# indices from the tri.simplices command (line 27), pertaining to the points array, as the same indices
# for the spins at the corresponding sites.

subwind = np.zeros(sublattices)
for i in range(sublattices):
    # Definition holds only for 2D systems, so you must select a layer "slice" to take of your systems.
    # Here, I select the plane z=0.
    if sublattices == 1:
        layer = np.where(allpoints[:,2] == 1*2.715)[0]
    elif sublattices > 1:
        layer = np.where(np.logical_and(allpoints[:,2] == 0, atominfo[:,0] == i))[0]
    points = allpoints[layer]
    spins = allspins[layer]
    # Here, I remove the z-axis since my plane is for constant z.
    points = np.delete(points, 2, axis=1)
    # Delauney triangulation to create a triangle mesh of all lattice sites.
    tri=Delaunay(points)
    wind = 0
    # Tri.simplices contains the information about what points (given by their index in the points array) belong
    # to each triangle group.
    row, col = tri.simplices.shape
    ptcharge = np.zeros(row)
    pts = np.zeros((row, 2))
    defect = np.zeros((row,2))
    for x in range (0,row):
        s1, s2, s3 = [tri.simplices[x,i] for i in range(0,3)]
        # Iterative summation of the terms for the discrete lattice winding number.
        pts[x] = (points[s1] + points[s2] + points[s3])/3
        ptcharge[x] = (2*np.arctan((np.dot(spins[s1], np.cross(spins[s2],spins[s3])))/(1 + np.dot(spins[s1], spins[s2]) + np.dot(spins[s2], spins[s3]) + np.dot(spins[s3], spins[s1]))))
        defect[x][0] = np.dot(spins[s1], np.cross(spins[s2], spins[s3])) # See discussion of presence of defects in system
        defect[x][1] = np.linalg.norm(spins[s1]+spins[s2]+spins[s3], ord = 2)
    subwind[i] = np.sum(ptcharge)/(4*np.pi)
    chargemin = np.min(ptcharge)
    chargemax = np.max(ptcharge)
    if np.abs(chargemin) < np.abs(chargemax):
        cbmax = np.abs(chargemin)
        cbmin = chargemin
    else:
        cbmax = chargemax
        cbmin = -1 * chargemax
    plt.scatter(pts[:,0], pts[:,1], c = ptcharge, cmap = 'coolwarm', s = 5, vmin = cbmin, vmax = cbmax)
    plt.triplot(points[:,0], points[:,1], tri.simplices, lw = 0.5, color = 'silver')
    plt.colorbar()
    plt.axis('equal')
    plt.axis('off')
    plt.savefig('charge.png')
    
    
wind = np.sum(subwind)
print(wind)
print(subwind)
