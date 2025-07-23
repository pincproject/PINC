#!/usr/bin/env python3
"""
PINC Object Creation Using HDF5 + VTK Export

Generates a sphere or box object mask in a 3D grid, writes it to HDF5,
then exports the same data as a VTK structured grid for ParaView.
"""
import h5py
import numpy as np
import os

import configparser

# PyVista for VTK export
try:
    import pyvista as pv
except ImportError:
    raise ImportError("PyVista is required for VTK export: pip install pyvista")

##################### INPUT ###########################
# PINC DOMAIN
domainSizeX = 32  # Domain size in X (cells)
domainSizeY = 32  # Domain size in Y (cells)
domainSizeZ = 32  # Domain size in Z (cells)

# PINC OBJECT
OBJ = "cubesat"     # "sphere", "box", or "cubesat"

# OBJECT LOCATION (center for sphere, origin for box)
ObjLocX = domainSizeX // 2
ObjLocY = domainSizeY // 2
ObjLocZ = domainSizeZ // 2

# CubeSat with booms parameters (used when OBJ="cubesat")
CubeSatSize = 4       # Cube side length in cells
BoomLength   = 6      # Boom length in cells (extends from cube face)
BoomThickness= 1       # Boom thickness in cells
# Boom directions: choose any of "+X", "-X", "+Y", "-Y", "+Z", "-Z"
# BoomDirs = ["+X", "-X", "+Y", "-Y", "+Z", "-Z"]
BoomDirs = ["+X", "-X"]

# Path to the debris.ini file (adjust as needed)
DebrisIni = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'input', 'debris.ini'))

if OBJ == "sphere":
    ObjR = 3        # Object radius (cells)
elif OBJ == "box":
    ObjSizeX = 4    # Box size in X (cells)
    ObjSizeY = 4    # Box size in Y (cells)
    ObjSizeZ = 4    # Box size in Z (cells)
elif OBJ == "cubesat":
    # Parameters defined above: CubeSatSize, BoomLength, BoomThickness
    pass
else:
    raise ValueError("OBJ must be 'sphere', 'box', or 'cubesat'")

# Create empty domain array: shape (Z, Y, X)
domain = np.zeros((domainSizeZ, domainSizeY, domainSizeX), dtype=int)

# Fill domain mask
if OBJ == "sphere":
    # Spherical object: sample in radial, azimuthal, polar coords
    rad, phi, theta = np.mgrid[
        0:ObjR:50j,
        0:2*np.pi:50j,
        0:np.pi:50j
    ]
    # Compute float coordinates
    Xf = ObjLocX + rad * np.sin(theta) * np.cos(phi)
    Yf = ObjLocY + rad * np.sin(theta) * np.sin(phi)
    Zf = ObjLocZ + rad * np.cos(theta)
    # Round, cast, and clamp
    X = np.clip(np.round(Xf).astype(int), 0, domainSizeX - 1)
    Y = np.clip(np.round(Yf).astype(int), 0, domainSizeY - 1)
    Z = np.clip(np.round(Zf).astype(int), 0, domainSizeZ - 1)
    # Vectorized fill
    domain[Z, Y, X] = 1

elif OBJ == "box":
    # Box object: fill axis-aligned block
    domain[
        ObjLocZ : ObjLocZ + ObjSizeZ,
        ObjLocY : ObjLocY + ObjSizeY,
        ObjLocX : ObjLocX + ObjSizeX
    ] = 1

elif OBJ == "cubesat":
    # Central cube
    half = CubeSatSize // 2
    z0 = max(0, ObjLocZ - half)
    z1 = min(domainSizeZ, ObjLocZ + half + CubeSatSize % 2)
    y0 = max(0, ObjLocY - half)
    y1 = min(domainSizeY, ObjLocY + half + CubeSatSize % 2)
    x0 = max(0, ObjLocX - half)
    x1 = min(domainSizeX, ObjLocX + half + CubeSatSize % 2)
    domain[z0:z1, y0:y1, x0:x1] = 1

    # Precompute boom cross‑section bounds
    yt0 = ObjLocY - BoomThickness//2
    yt1 = yt0 + BoomThickness
    zt0 = ObjLocZ - BoomThickness//2
    zt1 = zt0 + BoomThickness
    xt0 = ObjLocX - BoomThickness//2
    xt1 = xt0 + BoomThickness

    # Draw booms based on BoomDirs
    if "+X" in BoomDirs:
        bx0 = ObjLocX + half
        bx1 = min(domainSizeX, bx0 + BoomLength)
        domain[zt0:zt1, yt0:yt1, bx0:bx1] = 1
    if "-X" in BoomDirs:
        bx0 = max(0, ObjLocX - half - BoomLength)
        bx1 = ObjLocX - half
        domain[zt0:zt1, yt0:yt1, bx0:bx1] = 1
    if "+Y" in BoomDirs:
        by0 = ObjLocY + half
        by1 = min(domainSizeY, by0 + BoomLength)
        domain[zt0:zt1, by0:by1, xt0:xt1] = 1
    if "-Y" in BoomDirs:
        by0 = max(0, ObjLocY - half - BoomLength)
        by1 = ObjLocY - half
        domain[zt0:zt1, by0:by1, xt0:xt1] = 1
    if "+Z" in BoomDirs:
        bz0 = ObjLocZ + half
        bz1 = min(domainSizeZ, bz0 + BoomLength)
        domain[bz0:bz1, yt0:yt1, xt0:xt1] = 1
    if "-Z" in BoomDirs:
        bz0 = max(0, ObjLocZ - half - BoomLength)
        bz1 = ObjLocZ - half
        domain[bz0:bz1, yt0:yt1, xt0:xt1] = 1

# Ensure output directory exists
data_dir = 'data'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# 1) Write HDF5 file
h5_path = os.path.join(data_dir, 'object.grid.h5')
with h5py.File(h5_path, 'w') as hf:
    hf.create_dataset('Object', data=domain)
print(f"HDF5 written to: {h5_path}")

# 2) Write VTK file via PyVista
# Use ImageData since UniformGrid may not exist in this PyVista version
data = domain  # alias for clarity

# Build an ImageData (uniform grid) with point dimensions = (nx+1, ny+1, nz+1)
grid = pv.ImageData(
    dimensions=(domainSizeX + 1, domainSizeY + 1, domainSizeZ + 1),
    origin=(0.0, 0.0, 0.0),
    spacing=(1.0, 1.0, 1.0)
)
# Attach cell data (mask) flattened in C order
grid.cell_data['Object'] = data.flatten(order='C')

vtk_path = os.path.join(data_dir, 'object.vtk')
grid.save(vtk_path)
print(f"VTK written to: {vtk_path}")

# List contents of data directory
print("\nContents of data/:")
for fname in os.listdir(data_dir):
    fpath = os.path.join(data_dir, fname)
    size_kb = os.path.getsize(fpath) / 1024
    print(f"  {fname}\t{size_kb:.1f} KB")


# --- Compute surface area ---
# Read step size (in meters) from debris.ini
config = configparser.ConfigParser()
config.read(DebrisIni)
raw_step = config.get('grid', 'stepSize')
# take the first value before any semicolon
step_size_m = float(raw_step.split(';')[0])

# Cube surface area
cube_side = CubeSatSize * step_size_m
SA_cube = 6 * cube_side**2

# Booms surface area (4 side faces + 1 free end face per boom)
len_m = BoomLength * step_size_m
thick_m = BoomThickness * step_size_m
SA_booms = 0.0
for direction in BoomDirs:
    SA_booms += 4 * thick_m * len_m + thick_m**2

SA_total = SA_cube + SA_booms
print(f"\nStep size: {step_size_m:.6f} m")
print(f"Cube surface area: {SA_cube:.6f} m²")
print(f"Booms surface area: {SA_booms:.6f} m²")
print(f"Total spacecraft surface area: {SA_total:.6f} m²")