import numpy as np
import numpy.random as rd
from math import pi
from scipy.special import hyp2f1

HEIGHT_STEP  = 2.0   # use a grid with 2m resolution in z direction
GRID_STEP    = 2.0   # grid resolution in x and y direction
np.random.seed(1234) # do this only if you want to have reproducible wind fields

def pfq(z):
  return np.real(hyp2f1(1./3., 17./6., 4./3., z))

def createWindField(x, y, z,  sigma1 = None, gamma= 3.9, ae= 0.1, length_scale=33.6):
    """
    sigma1: Target value(s) for the turbulence
           intensity. This can either be a single value, in
           which case the statistics are corrected based on
           the longitudinal component only, or a vector where
           each u,v,w-component can be defined separately.
    gamma: wind shear (zero for isotropic turbulence, 3.9 for IEC wind profile)
    ae:    Coefficient of the inertial cascade, ae = a*e^(2/3),
           where a = 1.7 is the three-dimensional Kolmogorov
           constant and e = the mean TKE dissipation rate [m/s].
    Performance: Python:  17 seconds for a field of 50x800x200 m with 2 m resolution
                 Matlab:  17 seconds
    """
    if sigma1 is not None:
        if type(sigma1) != int and type(sigma1) != float and type(sigma1) != np.float64:
            if len(sigma1) != 3:
                raise Exception("The parameter 'sigma' must either be a single value or a 3-component vector.")
    if (x.ravel()[0] > x.ravel()[-1]) or (y.ravel()[0] > y.ravel()[-1]) or (z.ravel()[0] > z.ravel()[-1]):
        raise Exception("The values of x, y, and z must be monotonically increasing.")

    # Standard deviations
    sigma_iso = 0.55 * sigma1
    sigma2 = 0.7 * sigma1
    sigma3 = 0.5 * sigma1

    # Domain size
    # Number of elements in each direction
    nx = x.shape[0]  #size(x,2);
    ny = y.shape[1]  #size(y,1);
    nz = z.shape[2]  #size(z,3);

    X, Y, Z = x, y, z

    Lx = X.ravel()[-1] - X.ravel()[0]
    Ly = Y.ravel()[-1] - Y.ravel()[0]
    Lz = Z.ravel()[-1] - Z.ravel()[0]

    # Wave number discretization
    # Shifted integers
    y_range=np.linspace(-ny/2., ny/2. - 1,  num=ny)
    x_range=np.linspace(-nx/2., nx/2. - 1,  num=nx)
    z_range=np.linspace(-nz/2., nz/2. - 1,  num=nz)
    m2, m1, m3 = np.meshgrid(y_range, x_range, z_range)

    m1 = np.fft.ifftshift(m1 + 1e-6)
    m2 = np.fft.ifftshift(m2 + 1e-6)
    m3 = np.fft.ifftshift(m3 + 1e-6)

    # Wave number vectors
    k1 = 2 * pi * m1 * (length_scale / Lx)
    k2 = 2 * pi * m2 * (length_scale / Ly)
    k3 = 2 * pi * m3 * (length_scale / Lz)

    k = np.sqrt(k1**2 + k2**2 + k3**2)

    # Non-dimensional distortion time
    pfq_term = pfq(-k**-2)
    beta = gamma / (k**(2./3.) * np.sqrt(pfq_term))

    # Initial wave vectors
    k30 = k3 + beta * k1
    k0 = np.sqrt(k1**2 + k2**2 + k30**2)

    # Non-dimensional von Karman isotropic energy spectrum
    E0 = ae * length_scale**(5./3.) * k0**4 / (1 + k0**2)**(17./6.)

    # Correlation matrix
    C1 = (beta * k1**2 * (k1**2 + k2**2 - k3 * (k3 + beta * k1))) / (k**2 * (k1**2 + k2**2))
    C2 = (k2 * k0**2) / (k1**2 + k2**2)**(3./2.) * np.arctan2((beta * k1 * np.sqrt(k1**2 + k2**2)).real, (k0**2 \
         - (k3 + beta * k1) * k1 * beta).real)

    zeta1 = C1 - k2 / k1 * C2
    zeta2 = C2 + k2 / k1 * C1
    B = sigma_iso * np.sqrt(2.0 * pi**2 * length_scale**3 * E0 / (Lx * Ly *Lz * k0**4))
    # B = sigma_iso * sqrt(2*pi^2 * L^3 * E0 ./ (Lx*Ly*Lz * k0.^4));

    C = np.zeros((3, 3,) + X.shape)
    C[0,0] = B * k2 * zeta1
    C[0,1] = B * (k30 - k1 * zeta1)
    C[0,2] = B * -k2
    C[1,0] = B * (k2 * zeta2 - k30)
    C[1,1] = B * -k1 * zeta2
    C[1,2] = B * k1
    C[2,0] = B *  k2 * k0**2 / k**2
    C[2,1] = B * -k1 * k0**2 / k**2

    # Set up stochastic field
    # White noise vector
    n = rd.normal(0, 1, [3, 1] + list(k.shape)) + 1j * rd.normal(0, 1, [3, 1] + list(k.shape))

    # Stochastic field
    dZ = np.zeros(((3,) + k.shape), dtype=np.dtype(np.complex128))
    for i in range(X.shape[0]):
        for j in range(Y.shape[1]):
            for k in range(Z.shape[2]):
                dZ[:, i, j, k] = (np.dot(C[:, : , i, j, k] , n[:, :, i, j, k])).reshape((3,))

    # Reconstruct time series
    u = nx * ny * nz * (np.fft.ifftn(np.squeeze(dZ[0,:,:,:]))).real
    v = nx * ny * nz * (np.fft.ifftn(np.squeeze(dZ[1,:,:,:]))).real
    w = nx * ny * nz * (np.fft.ifftn(np.squeeze(dZ[2,:,:,:]))).real

    if sigma1 is not None:
        su = np.std(u)
        sv = np.std(v)
        sw = np.std(w)
        u = sigma1/su * u
        v = sigma2/sv * v
        w = sigma3/sw * w
    return u, v, w

def createGrid(ny=50, nx=100, nz=50, z_min=25, res=GRID_STEP):
    """
    res: resolution of the grid in x and y direction in meters
    ny:  number of meters in y direction
    nx:  number of meters in x direction (downwind)
    nz:  number of meters in z direction (up)
    z_min: minimal height in m
    """
    res = int(res)
    y_range=np.linspace(-ny//2, ny//2,      num=ny//res+1)
    x_range=np.linspace(0, nx,            num=nx//res+1)
    z_range=np.linspace(z_min, z_min+nz,  num=nz//int(HEIGHT_STEP)+1)
    # returns three 3-dimensional arrays with the components of the position of the grid points
    y, x, z = np.meshgrid(y_range, x_range, z_range)
    return y, x, z

y, x, z = createGrid(10, 20, 10, 5)
u, v, w = createWindField(x, y, z, sigma1=1)
print(u.shape)
print(v.shape)
print(w.shape)
