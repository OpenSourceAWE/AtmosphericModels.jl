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

