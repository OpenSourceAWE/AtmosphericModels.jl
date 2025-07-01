import numpy as np
from math import pi, sqrt

GRID_STEP   = 2.0    # Resolution of the grid in x and y direction in meters
HEIGHT_STEP = 2.0    # Resolution in z direction in meters

def create_grid(ny=50, nx=100, nz=50, z_min=25, res=GRID_STEP):
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

y, x, z = create_grid(100, 4050, 500, 70)

