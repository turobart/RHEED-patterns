import numpy as np
import matplotlib.pyplot as plt
import math
import sympy as sym
from sympy.utilities.lambdify import lambdify
import time

start_time = time.time()
# variables of the material
xSn = 0.3
# latticeConstant = (6.1240-xSn*0.1246) #in the given direction
latticeConstant = 6.490 #aSn
theta_deg = 0.
theta = np.deg2rad(theta_deg)

# broadening of the streaks to produce realistic RHEED image
sigma_x=0.1
sigma_y=0.1/3*2

# variables of the system
d = 32 #distance between sample and RHEED screen in cm. 32 cm in GM#2 VEECO GENxplor MBE System 
incAngle = np.deg2rad(3.2) # angle of incidence of the electron beam. Typical range: 2-5 degree
acceleratingV = 15000 #accelerating voltage of the electron beam source. Typiucal range 10,000-15,000 eV

electronWavelength = 12.247/math.sqrt(acceleratingV*(1+acceleratingV*0.000001))
k0=2*math.pi/electronWavelength # radius of the Ewald sphere

sampling_size = 500 # 500 is suggested for initial imageing; more than 1500 is suggested for fine images

# size of G reciprocal array
N=10.
M=6.
# creation of the array of atoms
nv=(np.arange(-N,N+1))
mv=(np.arange(-N,M+1))

radius_shift = math.cos(incAngle)*k0 # shift of the Ewald sphere in the reciprocal space

RotMatrix = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
# boudaries of the k-space calculations. 
# calculated using 
# x = k0*4//math.sqrt(d**2+4**2); y = k0*(-d/math.sqrt(d**2+5**2)+math.cos(incAngle))
# 4 and 5 correspond to the dimensions of the real RHEED screen
kxmin, kxmax =  -8, 8
kymin, kymax =  radius_shift-k0, 1.25

def defineLattice(latticeA1, latticeA2, latticeUnitCellFactor):
    latticeA1 = (latticeA1)*latticeConstant*latticeUnitCellFactor
    latticeA2 = (latticeA2)*latticeConstant*latticeUnitCellFactor
    
    # calculation of the reciprocal space distances
    latticeB1 = ((2*math.pi/np.cross(latticeA1, latticeA2))*np.array([latticeA2[1], -latticeA2[0]])).round(5)
    latticeB2 = ((2*math.pi/np.cross(latticeA1, latticeA2))*np.array([-latticeA1[1], latticeA1[0]])).round(5)
    
    # creation of the reciprocal space array
    latticeGxp=np.array([[i*latticeB1[0]+j*latticeB2[0] for i in nv] for j in mv]).round(4)
    latticeGyp=np.array([[i*latticeB1[1]+j*latticeB2[1] for i in nv] for j in mv]).round(4)

    electron_scattering_factor=[1,1] # array of the electron scattering factors. Must have lengh of aM.
    atomPosition=[0,1] # array of the relative atoms positions in the unit cell. Must have lengh of aM.
    aM=2
    
    # generation of grid for image calculation
    kx = np.linspace(kxmin, kxmax, sampling_size)
    ky = np.linspace(kymin, kymax, sampling_size)
    xx, yy = np.meshgrid(kx, ky)
    # rotation
    xxRot, yyRot = np.einsum('ij, mni -> jmn', RotMatrix, np.dstack([xx, yy]))
    
    # generation of the reciprocal lattice rods projected on 2D surface
    RLRP = np.zeros((sampling_size, sampling_size))
    for h in range(len(mv)):
        for k in range(len(nv)):       
            structureFactor = sum((electron_scattering_factor[aM-1]*np.exp(((latticeGxp[h][k]+latticeGyp[h][k])*atomPosition[aM-1])*1j)) for i in range(1, aM + 1))
            DiffIntensity = np.square(np.absolute(structureFactor))
            RLRP += DiffIntensity*np.exp(-((xxRot-latticeGxp[h][k])**2/(sigma_x**2) + (yyRot-latticeGyp[h][k])**2/(sigma_y**2)))
            
    return RLRP
    
def calculateRHEEDarray(RLRP):
    x,y, xe, ye = sym.symbols('x,y, xe, ye', positive = True)
    dd, kk, aa = sym.symbols('dd kk aa', positive = True)
       
    kx_c = np.linspace(kxmin, kxmax, sampling_size, dtype=complex)
    ky_c = np.linspace(kymin, kymax, sampling_size, dtype=complex)
    xx, yy = (np.meshgrid(kx_c, ky_c))

#    # solution calculation with sympy
#    # ----------
#     Eq1 = sym.Eq(kk*x/sym.sqrt(dd**2+x**2+y**2), xe)
#     Eq2 = sym.Eq(kk*(-dd/sym.sqrt(dd**2+x**2+y**2)+sym.cos(aa)), ye)
#     sol = sym.solve([Eq1,Eq2],(x,y))
#     print(sol)
    # ----------
    eqs = (-sym.sqrt(dd**2*xe**2)/(kk*sym.cos(aa) - ye), dd*sym.sqrt(kk**2*sym.sin(aa)**2 + 2*kk*ye*sym.cos(aa) - xe**2 - ye**2)/(kk*sym.cos(aa) - ye))
    values = {kk: k0, dd: d, aa: incAngle}
    solsol = sym.Tuple(eqs).subs(values).simplify()
    
    func = lambdify((xe,ye), solsol,'numpy') # numpy-ready function for array calculations

    numpy_array_of_results = np.real((func(xx, yy)))

    # generate RHEED image
    plot_name = 'Angle_%.2f' %theta_deg
    plot_numpy = plt.figure(plot_name)
    plt.axes().set_aspect('equal')
    plt.pcolormesh(numpy_array_of_results[0][0], numpy_array_of_results[0][1], RLRP, cmap='gray')
    plt.ylim(0, 5) # in cm, should correspond to real size of RHEED screen 
    plt.xlim(-4, 4) # in cm, should correspond to real size of RHEED screen
    # image saving
    if theta_deg%1>0:
        fig_name = 'Angle_%.2f' %theta_deg
        fig_name += '.png'
        print(fig_name)
    else:
        fig_name = 'Angle_%d' %theta_deg +'.png'
    plot_numpy.savefig(fig_name, bbox_inches='tight',transparent=True,format='png')

# definition of the basis vectors; uncomment as needed
#growth direction (111)
a11 = np.array([math.sqrt(3)/2, -1/2])
a12 = np.array([0., 1.])
# a21 = np.array([math.sqrt(3)/2, -1/2])
# a22 = np.array([0., 1.])
#growth direction (001)
# a11 = np.array([2., 0.])
# a12 = np.array([0., 1.])
# a21 = np.array([1., 0.])
# a22 = np.array([0., 2.])

# calculation of atoms distances in topmost layer
UnitCellFactor = math.sqrt(2)/2 # used to take into account the disribution of atoms on the topmost layer.
                                # Derived from geometrical properties of the unit cell and desired growt orientation

RLRP1 = defineLattice(a11, a12, UnitCellFactor)
# IDF2 = defineLattice(a21, a22, UnitCellFactor)

RLRP = RLRP1# + RLRP2
        
RLRP = RLRP/np.max(RLRP)
RLRP = np.array(RLRP * 255, dtype = np.uint8)

calculateRHEEDarray(RLRP)
print("--- %s seconds ---" % (time.time() - start_time))
plt.show()