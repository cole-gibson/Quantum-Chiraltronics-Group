This is a utility to calculate the magnetocrystalline anisotropy coefficients for a magnetic material from calculations of the total energy. There are two requisite files:

  1) input: The first line must contain the Schoenflies symbol for the crystal point group of the material. The following lines contain the crystal basis, with each vector wrapped as [...] and belonging to a new line. These should be presented in the order a, b, and c based on the conventional choices of these vectors.
  2) energies: This file must contain the theta, phi, and energy values. Here, theta and phi are the azimuthal and polar angles of the magnetization the energy was calculated for respectively.

The fitting is done using expansions derived in this paper: Döring, W. Die Richtungsabhängigkeit der Kristallenergie. Ann. Phys. 456, 102–109 (1957).

The output includes the expansion used, easy axis, hard axis, maximum energy value, and an equation in spherical coordinates for the energy surface with lowest energy = 0.
