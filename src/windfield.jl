"""
This module provides methods to create a turbulent wind field and to read 
the actual wind velocity vector as function of the 3d position
vector and of the time.
It additionally provides functions for plotting the rel. turbulence as function
of the height.
This code is based on the Matlab module Mann.m written by René Bos.
The code is based on the following papers:
 - Mann, Jakob (1994). The Spatial Structure of Neutral Atmospheric
   Surface-Layer Turbulence. Journal of Fluid Mechanics 273,
   pp. 141-168.
 - Mann, Jakob. (1998). Wind Field Simulation. Probabilistic
   Engineering Mechanics 13(4), pp. 269-282.
"""

function pfq(z)
    _₂F₁(1. /3 , 17. /6, 4. /3, z)
end

function calc_sigma1(am, v_wind_gnd)
    v_height = v_wind_gnd * calc_wind_factor(am, am.set.avg_height, Val{Int(EXP)})
    am.set.i_ref * (0.75 * v_height + 5.6)
end

"""
Find 2^n that is equal to or greater than i.
"""
function nextpow2(i)
    n = 1
    while n < i
        n *= 2
    end
    n
end

using NPZ

const HOME = ENV["HOME"]  # Or set this manually
const V_WIND_GND = 8.0    # Default value, change as needed

function calcFullName(v_wind_gnd, basename="windfield_4050_500", rel_sigma=1.0)
    path = get_data_path() * "/"
    name = basename * "_" * @sprintf("%.1f", rel_sigma)
    name *= "_" * @sprintf("%.1f", v_wind_gnd)
    return path * name
end

function save(x, y, z, u, v, w, param; basename="windfield_4050_500", v_wind_gnd=V_WIND_GND)
    fullname = calcFullName(v_wind_gnd, basename=basename)
    # Save as compressed .npz
    NPZ.npzwrite(fullname * ".npz", Dict(
        "x" => x,
        "y" => y,
        "z" => z,
        "u" => u,
        "v" => v,
        "w" => w,
        "param" => param
    ))
end

function load(; basename="windfield_4050_500", v_wind_gnd=8.0)
    fullname = calcFullName(v_wind_gnd, basename=basename)
    npzfile = NPZ.npzread(fullname * ".npz")
    return npzfile["x"], npzfile["y"], npzfile["z"], npzfile["u"], npzfile["v"], npzfile["w"], npzfile["param"]
end

# def loadWindField(speed):
#     """
#     1. determine the nearest wind speed from the list V_WIND_GNDS
#     2. load the wind field for this wind speed
#     """
#     if abs(speed - V_WIND_GNDS[1]) < abs(speed - V_WIND_GNDS[2]) and abs(speed - V_WIND_GNDS[1]) < abs(speed - V_WIND_GNDS[0]):
#     # if abs(speed - 5.2) < abs(speed - 8.0) and abs(speed - 5.2) < abs(speed - 3.8):
#         return load(v_wind_gnd = V_WIND_GNDS[1])
#     elif abs(speed - V_WIND_GNDS[2]) < abs(speed - V_WIND_GNDS[1]) and abs(speed - V_WIND_GNDS[2]) < abs(speed - V_WIND_GNDS[0]):
#         return load(v_wind_gnd = V_WIND_GNDS[2])
#     else:
#         return load(v_wind_gnd = V_WIND_GNDS[0])

# def createGrid(ny=50, nx=100, nz=50, z_min=25, res=GRID_STEP):
#     """
#     res: resolution of the grid in x and y direction in meters
#     ny:  number of meters in y direction
#     nx:  number of meters in x direction (downwind)
#     nz:  number of meters in z direction (up)
#     z_min: minimal height in m
#     """
#     y_range=np.linspace(-ny/2, ny/2,      num=ny/res+1)
#     x_range=np.linspace(0, nx,            num=nx/res+1)
#     z_range=np.linspace(z_min, z_min+nz,  num=nz/HEIGHT_STEP+1)
#     # returns three 3-dimensional arrays with the components of the position of the grid points
#     y, x, z = np.meshgrid(y_range, x_range, z_range)
#     return y, x, z

# def showGrid(x, y, z):
#     """
#     x: downwind direction
#     z: up
#     y: orthogonal to x and z
#     """
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.scatter(x, y, z)

#     ax.set_xlabel('X Label')
#     ax.set_ylabel('Y Label')
#     ax.set_zlabel('Height [m]')

# def show2Dfield(X, Y, U, V, scale=1.0, width=0.07):
#     """
#     x, y: position of the wind velocity vectors
#     u, v: wind velocity vector components
#     """
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.quiver(X, Y, U, V, angles='uv', units='xy', scale=scale, width=width)

# def plotTurbulenceVsHeight(X, Y, Z, U, V, W):
#     fig = plt.figure()
#     height_min = np.min(Z)
#     height_max = np.max(Z)
#     print "height range: ", height_min, height_max
#     Height, XX, YY, ZZ = [], [], [], []
#     for height in np.linspace(height_min, height_max, int(round((height_max-height_min) / HEIGHT_STEP))+1):
#         Height.append(height)
#         v_wind = calcWindHeight(V_WIND_GND, height)
#         height_index = int(round((height-height_min) / HEIGHT_STEP))
#         u2 =  U[:,:,height_index]; XX.append(u2.std() / v_wind * 100.0)
#         v2 =  V[:,:,height_index]; YY.append(v2.std() / v_wind * 100.0)
#         w2 =  W[:,:,height_index]; ZZ.append(w2.std() / v_wind * 100.0)
#     plt.plot(Height, XX, label = 'rel. turbulence x [%]', color="black")
#     plt.plot(Height, YY, label = 'rel. turbulence y [%]', color="blue")
#     plt.plot(Height, ZZ, label = 'rel. turbulence z [%]', color="green")
#     plt.legend(loc='upper right')
#     plt.gca().set_ylabel('Rel. turbulence [%]')
#     plt.gca().set_xlabel('Height [m]')

# def createWindField(x, y, z,  sigma1 = None, gamma= 3.9, ae= 0.1, length_scale=33.6):
#     """
#     sigma1: Target value(s) for the turbulence
#            intensity. This can either be a single value, in
#            which case the statistics are corrected based on
#            the longitudinal component only, or a vector where
#            each u,v,w-component can be defined separately.
#     gamma: wind shear (zero for isotropic turbulence, 3.9 for IEC wind profile)
#     ae:    Coefficient of the inertial cascade, ae = a*e^(2/3),
#            where a = 1.7 is the three-dimensional Kolmogorov
#            constant and e = the mean TKE dissipation rate [m/s].
#     Performance: Python:  17 seconds for a field of 50x800x200 m with 2 m resolution
#                  Matlab:  17 seconds
#     """
#     if sigma1 is not None:
#         if type(sigma1) != int and type(sigma1) != float and type(sigma1) != np.float64:
#             if len(sigma1) != 3:
#                 raise Exception("The parameter 'sigma' must either be a single value or a 3-component vector.")
#     if (x.ravel()[0] > x.ravel()[-1]) or (y.ravel()[0] > y.ravel()[-1]) or (z.ravel()[0] > z.ravel()[-1]):
#         raise Exception("The values of x, y, and z must be monotonically increasing.")

#     # Standard deviations
#     sigma_iso = 0.55 * sigma1
#     sigma2 = 0.7 * sigma1
#     sigma3 = 0.5 * sigma1

#     # Domain size
#     # Number of elements in each direction
#     nx = x.shape[0]  #size(x,2);
#     ny = y.shape[1]  #size(y,1);
#     nz = z.shape[2]  #size(z,3);

#     X, Y, Z = x, y, z

#     Lx = X.ravel()[-1] - X.ravel()[0]
#     Ly = Y.ravel()[-1] - Y.ravel()[0]
#     Lz = Z.ravel()[-1] - Z.ravel()[0]

#     # Wave number discretization
#     # Shifted integers
#     y_range=np.linspace(-ny/2., ny/2. - 1,  num=ny)
#     x_range=np.linspace(-nx/2., nx/2. - 1,  num=nx)
#     z_range=np.linspace(-nz/2., nz/2. - 1,  num=nz)
#     m2, m1, m3 = np.meshgrid(y_range, x_range, z_range)

#     m1 = np.fft.ifftshift(m1 + 1e-6)
#     m2 = np.fft.ifftshift(m2 + 1e-6)
#     m3 = np.fft.ifftshift(m3 + 1e-6)

#     # Wave number vectors
#     k1 = 2 * pi * m1 * (length_scale / Lx)
#     k2 = 2 * pi * m2 * (length_scale / Ly)
#     k3 = 2 * pi * m3 * (length_scale / Lz)
#     k = np.sqrt(k1**2 + k2**2 + k3**2)

#     # Non-dimensional distortion time
#     pfq_term = pfq(-k**-2)
#     beta = gamma / (k**(2./3.) * np.sqrt(pfq_term))

#     # Initial wave vectors
#     k30 = k3 + beta * k1
#     k0 = np.sqrt(k1**2 + k2**2 + k30**2)

#     # Non-dimensional von Karman isotropic energy spectrum
#     E0 = ae * length_scale**(5./3.) * k0**4 / (1 + k0**2)**(17./6.)

#     # Correlation matrix
#     C1 = (beta * k1**2 * (k1**2 + k2**2 - k3 * (k3 + beta * k1))) / (k**2 * (k1**2 + k2**2))
#     C2 = (k2 * k0**2) / (k1**2 + k2**2)**(3./2.) * np.arctan2((beta * k1 * np.sqrt(k1**2 + k2**2)).real, (k0**2 \
#          - (k3 + beta * k1) * k1 * beta).real)

#     zeta1 = C1 - k2 / k1 * C2
#     zeta2 = C2 + k2 / k1 * C1
#     B = sigma_iso * np.sqrt(2.0 * pi**2 * length_scale**3 * E0 / (Lx * Ly *Lz * k0**4))
#     # B = sigma_iso * sqrt(2*pi^2 * L^3 * E0 ./ (Lx*Ly*Lz * k0.^4));

#     C = np.zeros((3, 3,) + X.shape)
#     C[0,0] = B * k2 * zeta1
#     C[0,1] = B * (k30 - k1 * zeta1)
#     C[0,2] = B * -k2
#     C[1,0] = B * (k2 * zeta2 - k30)
#     C[1,1] = B * -k1 * zeta2
#     C[1,2] = B * k1
#     C[2,0] = B *  k2 * k0**2 / k**2
#     C[2,1] = B * -k1 * k0**2 / k**2

#     # Set up stochastic field
#     # White noise vector
#     n = rd.normal(0, 1, [3, 1] + list(k.shape)) + 1j * rd.normal(0, 1, [3, 1] + list(k.shape))

#     # Stochastic field
#     dZ = np.zeros(((3,) + k.shape), dtype=np.dtype(np.complex128))
#     for i in range(X.shape[0]):
#         for j in range(Y.shape[1]):
#             for k in range(Z.shape[2]):
#                 dZ[:, i, j, k] = (np.dot(C[:, : , i, j, k] , n[:, :, i, j, k])).reshape((3,))

#     # Reconstruct time series
#     u = nx * ny * nz * (np.fft.ifftn(np.squeeze(dZ[0,:,:,:]))).real
#     v = nx * ny * nz * (np.fft.ifftn(np.squeeze(dZ[1,:,:,:]))).real
#     w = nx * ny * nz * (np.fft.ifftn(np.squeeze(dZ[2,:,:,:]))).real

#     if sigma1 is not None:
#         su = np.std(u)
#         sv = np.std(v)
#         sw = np.std(w)
#         u = sigma1/su * u
#         v = sigma2/sv * v
#         w = sigma3/sw * w
#     return u, v, w

# def addWindSpeed(z, u, v_wind_ground=V_WIND_GND):
#     """
#     Modify the velocity component u such that the average wind speed, calculated according
#     to the given wind profile is added.
#     """
#     min_height = np.min(z)
#     max_height = np.max(z)
#     print "Minimal height: ", min_height
#     print "Maximal height: ", max_height
#     for i in range(z.shape[2]):
#         height = z[0,0,i]
#         v_wind = calcWindHeight(v_wind_ground, height)
#         # print i, v_wind
#         u[:,:,i] += v_wind

# class WindField(object):
#     def __init__(self, speed = V_WIND_GND):
#         try:
#             self.last_speed = 0.0
#             self.load(speed)
#             self.valid = True
#         except:
#            self.valid = False
#            print "Error reading wind field!"

#     def load(self, speed):
#         global ALPHA, V_WIND_GND
#         if speed == self.last_speed:
#             return
#         print "Loading wind field ...", form(speed)
#         self.last_speed = speed
#         self.x, self.y, self.z, self.u, self.v, self.w, param = loadWindField(speed)
#         self.x_max = np.max(self.x)
#         self.y_max = np.max(self.y)
#         self.x_min = np.min(self.x)
#         self.y_min = np.min(self.y)
#         self.y_range = self.y_max - self.y_min
#         self.x_range = self.x_max - self.x_min
#         self.u_pre, self.v_pre, self.w_pre = [ndimage.spline_filter(item, order=3) for item \
#                                                                                    in [self.u, self.v, self.w]]
#         ALPHA = param[0]
#         V_WIND_GND = param[1]
#         print "Finished loading wind field."

#     def getVWindGnd(self):
#         return V_WIND_GND

#     def setVWindGnd(self, v_wind_gnd):
#         global V_WIND_GND
#         V_WIND_GND = v_wind_gnd
#         return

#     def getWind(self, x, y, z, t, interpolate=True, rel_turb = 0.351):
#         """ Return the wind vector for a given position and time. Linear interpolation in x, y and z.
#         3.34 sec for order = 2; 37µs for order = 1
#         """
#         assert(z >= 5.)
#         if z < 10.:
#             z = 10.0
#         assert(t >= 0.)
#         while x < 0.0:
#             x += self.x_range
#         while y > self.y_max:
#             y -= self.y_range
#         while y < self.y_min:
#             y += self.y_range
#         y1 = ((y + self.y_range / 2) / GRID_STEP)
#         v_wind_height = calcWindHeight(V_WIND_GND, z)
#         # print "--->", v_wind_height
#         # sys.exit()
#         x1 = (x + t * v_wind_height) / GRID_STEP
#         while x1 > self.u.shape[0] - 1:
#             x1 -= self.u.shape[0] - 1
#         z1 = z / HEIGHT_STEP
#         if z1 > self.u.shape[2] - 1:
#             z1 = self.u.shape[2] - 1
#         if z1 < 0:
#             z1 = 0
#         if interpolate:
#             x_wind = ndimage.map_coordinates(self.u_pre, [[x1], [y1], [z1]], order=3, prefilter=False)
#             y_wind = ndimage.map_coordinates(self.v_pre, [[x1], [y1], [z1]], order=3, prefilter=False)
#             z_wind = ndimage.map_coordinates(self.w_pre, [[x1], [y1], [z1]], order=3, prefilter=False)
#             v_x, v_y, v_z = x_wind[0]*rel_turb + v_wind_height, y_wind[0]*rel_turb, z_wind[0]*rel_turb
#             # v_x, v_y, v_z = v_wind_height, 0.0, 0.0
#             v_wind = sqrt(v_x*v_x + v_y*v_y + v_z*v_z)
#             if v_wind < 1.0:
#                 print "x1, y1, z1", x1, y1, z1
#             # print " v_x, v_y, v_z",  v_x, v_y, v_z
#             return v_x, v_y, v_z
#         else:
#             vx, vy, vz = self.u[x1, y1, z1] + v_wind_height, self.v[x1, y1, z1], self.w[x1, y1, z1]
#             return vx, vy, vz

# WIND_FIELD = WindField()

# def plotWindVsTime(x=0.0, y=0.0, z=197.3, rel_turb = REL_TURB[TEST]):
#     print "Relative turbulence: ", rel_turb
#     fig = plt.figure()
#     TIME, v_wind_x, v_wind_norm = [], [], []
#     for t in np.linspace(0.0, 600.0, 600*20):
#         TIME.append(t)
#         v_x, v_y, v_z = WIND_FIELD.getWind(x, y, z, t, rel_turb = rel_turb)
#         # print " v_x, v_y, v_z",  v_x, v_y, v_z
#         # sys.exit()
#         v_wind = sqrt(v_x*v_x + v_y*v_y + v_z*v_z)
#         v_wind_norm.append(v_wind)
#         v_wind_x.append(v_x)
#         if v_wind < 0.1:
#             print "Error for x, y, z, t: ", x, y, z, t
#         # print t, v_x
#     v_wind_x = np.array(v_wind_x)
#     su = np.std(v_wind_x)
#     mean =  v_wind_x.mean()
#     print "Mean wind x, standart deviation, turbulence intensity [%]: ", form(mean), form(su), \
#                                                                          form(su/mean * 100.0)
#     plt.plot(TIME, v_wind_x, label = 'Abs. wind speed at 197.3 m [m/s]', color="black")
#     # plt.legend(loc='upper right')
#     plt.gca().set_xlabel('Time [s]')
#     plt.gca().set_ylabel('Abs. wind speed at 197.3 m height [m/s]')

# def plotWindVsY(X, Y, Z, U, V, W, x=200, z=200.0):
#     fig = plt.figure()
#     y_min = np.min(Y)
#     y_max = np.max(Y)
#     print "y range: ", y_min, y_max
#     YY, TURB_1, TURB_2, TURB_3 = [], [], [], []
#     if False:
#         for y in np.linspace(y_min, y_max, int(round((y_max-y_min) / GRID_STEP)) + 1):
#             YY.append(y)
#             u2_1 = U[x, y, (z-100.) / HEIGHT_STEP]; TURB_1.append(u2_1)
#             u2_2 = U[x, y, z        / HEIGHT_STEP]; TURB_2.append(u2_2)
#             u2_3 = U[x, y, (z+100.) / HEIGHT_STEP]; TURB_3.append(u2_3)
#     else:
#         t = 0.0
#         for y in np.linspace(2*y_min, 2*y_max, 1000):
#             YY.append(y)
#             u2_1, v, w = WIND_FIELD.getWind(x, y, z-100., t) #U[x, y, (z-100.) / HEIGHT_STEP];
#             TURB_1.append(u2_1)
#             u2_2, v, w= WIND_FIELD.getWind(x, y, z, t) #U[x, y, z        / HEIGHT_STEP];
#             TURB_2.append(u2_2)
#             u2_3, v, w = WIND_FIELD.getWind(x, y, z+100., t) #U[x, y, (z+100.) / HEIGHT_STEP];
#             TURB_3.append(u2_3)
#     #YY.extend((y_max - y_min + GRID_STEP) + np.array(YY))
#     #TURB_1.extend(TURB_1); TURB_2.extend(TURB_2); TURB_3.extend(TURB_3)
#     plt.plot(YY, TURB_1, label = 'abs. wind x at 100m [m/s]', color="black")
#     plt.plot(YY, TURB_2, label = 'abs. wind x at 200m [m/s]', color="blue")
#     plt.plot(YY, TURB_3, label = 'abs. wind x at 300m [m/s]', color="red")
#     plt.legend(loc='upper right')
#     plt.gca().set_xlabel('Y position [m]')

# def newWindField(v_wind_gnd):
#     """
#     Create and save a new wind field for the given ground wind speed.
#     """
#     print "Creating wind field. This might take 10 minutes or more..."
#     y, x, z = createGrid(100, 4050, 500, 70)
#     # y, x, z = createGrid(50, 16200, 200, 100) # 600s at 27 m/s
#     # y, x, z = createGrid(10, 20, 10, 5)
#     if True:
#         sigma1 = REL_SIGMA * calcSigma1(v_wind_gnd)
#         u, v, w = createWindField(x, y, z, sigma1=sigma1)
#     else:
#         u, v, w = createWindField(x, y, z)
#     param = np.array((ALPHA, v_wind_gnd))
#     # addWindSpeed(z, u)
#     save(x, y, z, u, v, w, param, v_wind_gnd=v_wind_gnd)
#     print "Finshed creating and saving wind field!"

# def newWindFields():
#     """
#     Create and save new wind fields for all ground wind speeds, defined in V_WIND_GNDS.
#     """

#     for v_wind_gnd in V_WIND_GNDS:
#         print "Creating wind field. This might take 10 minutes or more..."
#         y, x, z = createGrid(100, 4050, 500, 70)
#         # y, x, z = createGrid(50, 16200, 200, 100) # 600s at 27 m/s
#         # y, x, z = createGrid(10, 20, 10, 5)
#         if True:
#             sigma1 = REL_SIGMA * calcSigma1(v_wind_gnd)
#             u, v, w = createWindField(x, y, z, sigma1=sigma1)
#         else:
#             u, v, w = createWindField(x, y, z)
#         param = np.array((ALPHA, v_wind_gnd))
#         # addWindSpeed(z, u)
#         save(x, y, z, u, v, w, param, v_wind_gnd=v_wind_gnd)
#         print "Finshed creating and saving wind field!"
#         del y, x, z, u, v, w

# if __name__ == "__main__":
#     SAVE = False # True: calculate and save new wind field; False: use saved wind field
#     if not SAVE:
#         if WIND_FIELD.valid:
#             x, y, z = WIND_FIELD.x, WIND_FIELD.y, WIND_FIELD.z
#             u, v, w = WIND_FIELD.u, WIND_FIELD.v, WIND_FIELD.w
#         else:
#             SAVE = True
#     v_wind = calcWindHeight(V_WIND_GND, 200)
#     if SAVE:
#         newWindField(V_WIND_GND)
#         WIND_FIELD = WindField()

#     if False:
#         showGrid(x, y, z)
#     if False:
#         plotTurbulenceVsHeight(x, y, z, u, v, w)
#     if False:
#         plotWindVsY(x, y, z, u, v, w)
#     if False:
#         u2 =  u[:,1,:]
#         w2 =  w[:,1,:]
#         show2Dfield(x[:,1,:], z[:,1,:], u2, w2, scale=8.0)
#         # showGrid(X, Y, Z)
#     if True:
#         v_x, v_y, v_z = WIND_FIELD.getWind(0, 0, 197, 0)
#         print v_x, v_y, v_z
#         plotWindVsTime(0., 0., 197.)
#     if True:
#         print "sigma1", REL_TURB[2] * calcSigma1(V_WIND_GNDS[TEST])
#     if True:
#         print "V_WIND_GND at 6 m", V_WIND_GNDS[TEST]
#         print "wind at 197m:", calcWindHeight(V_WIND_GNDS[TEST], 197.0)

