# RHEED-patterns
This script uses geometrical transformation proposed by C. Liu et al.,[2022](https://doi.org/10.1116/6.0001899) adapted to the Python language.
The code can be used to generate images of RHEED patterns at any azimuth angle using the user-defined latice constant of the material, the azimuthal rotation of the sample, the distance between the sample and the RHEED screen (d), the angle of incidence of the electron beam (θ) and the accelerating voltage of the electron beam source.
---
Structure of the script:
1. An array of atoms is generated using a user-defined set of basis vectors.
2. The reciprocal space array is calculated using array from point 1.
   Each point of the calculated array corresponds to a top view of an infinite reciprocal space rod perpendicular to the surface of the material.
3. Azimuthal rotation of the array from point 2 is performed.
4. A 2D Gaussian function (with user-defined broadening) is used at each point of the reciprocal space array to simulate finite thickness of the rods.
   The intensity of the points depends on the structure factor of the atoms.The resulting 2D image is a sum of the Gaussians at positions of reciprocal space rods.
5. The 2D RLRP image is transformed using set of equations:
   $$x = {{k_0 x_d } \over {\sqrt{d^2+{x_d}^2+{y_d}^2}}}$$
   $$y = {k_0(-{d} \over {\sqrt{d^2+{x_d}^2+{y_d}^2}+cosθ)}}$$
   where x and y are coordinates of the reciprocal space of the sample surface, k0 is the wave vector of the electrons, xd and yd are the coordinates of the RHEED image on the fluorescent screen
7. The simulated RHEED pattern image is displayed and saved using azimuth angle as a filename.
