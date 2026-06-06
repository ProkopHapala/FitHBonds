import numpy as np
import matplotlib.pyplot as plt
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Plot H2O sp3 lone pair densities')
parser.add_argument('--cs', type=float, default=0.5,
                    help='s orbital coefficient (0 to 1). Default 0.5 (standard sp3). '
                         'Lower = more p character, more directional, density farther from nucleus.')
parser.add_argument('--zeta', type=float, default=2.2,
                    help='Slater orbital exponent. Default 2.2. Lower = more diffuse.')
parser.add_argument('--clamp', action='store_true',
                    help='Clamp negative amplitudes to zero before squaring (max(0,psi)^2)')
parser.add_argument('--power', type=float, default=1.0,
                    help='Angular power n for (cos)^n. Higher = sharper angular function. Default 1.0.')
parser.add_argument('--amplitude', action='store_true',
                    help='Plot raw amplitude (psi) instead of density (psi^2). Shows +/- regions.')
parser.add_argument('-o', '--output', type=str, default='H2O_sp3',
                    help='Output file prefix. Default: H2O_sp3')
args = parser.parse_args()

c_s = args.cs
zeta = args.zeta
clamp = args.clamp
power = args.power
amplitude = args.amplitude
output_prefix = args.output

print(f"Parameters: c_s={c_s:.3f}, c_p={np.sqrt(1-c_s**2):.3f}, zeta={zeta:.2f}, clamp={clamp}, power={power:.1f}, amplitude={amplitude}")

# H2O geometry (from xyz file)
O  = np.array([ 0.0,          0.0,          0.0        ])
H1 = np.array([-0.7580632005, 0.6358101311, 0.0        ])
H2 = np.array([ 0.7580632005, 0.6358101311, 0.0        ])

# Local coordinate system on Oxygen
# z-axis: perpendicular to molecular plane (normalized cross product of OH vectors)
v1 = H1 - O
v2 = H2 - O
z_axis = np.cross(v1, v2)
z_axis = z_axis / np.linalg.norm(z_axis)

# x-axis: bisector of H-O-H angle (in molecular plane)
x_axis = (v1 / np.linalg.norm(v1) + v2 / np.linalg.norm(v2))
x_axis = x_axis / np.linalg.norm(x_axis)

# y-axis: perpendicular to x and z (in molecular plane)
y_axis = np.cross(z_axis, x_axis)
y_axis = y_axis / np.linalg.norm(y_axis)

print(f"x-axis: {x_axis}")
print(f"y-axis: {y_axis}")
print(f"z-axis: {z_axis}")

# Helper: project global point to local frame
def to_local(r_global):
    dr = r_global - O
    return np.array([np.dot(dr, x_axis), np.dot(dr, y_axis), np.dot(dr, z_axis)])

# Slater-type orbitals (simplified, for visualization)
def orbital_s(r,  zeta=2.0):
    """1s-like radial function"""
    return np.exp(-zeta * r)

def orbital_px(r, x, y, z, zeta=2.0):
    return x * np.exp(-zeta * r)

def orbital_py(r, x, y, z, zeta=2.0):
    return y * np.exp(-zeta * r)

def orbital_pz(r, x, y, z, zeta=2.0):
    return z * np.exp(-zeta * r)

# Proper tetrahedral sp3 for water with C2v symmetry
# Local frame: x=bisector, y=in-plane perp, z=out-of-plane
# H atoms at ±52.25° from x-axis in xy plane (H-O-H angle = 104.5°)
# Lone pairs at ±52.25° from -x-axis in xz plane (pointing back, ±z)

# Tetrahedral half-angle: cos(52.25°) ≈ 0.612, sin(52.25°) ≈ 0.791
# These give 109.47° between any two sp3 directions

c = np.cos(np.radians(52.25))  # ≈ 0.612
s = np.sin(np.radians(52.25))  # ≈ 0.791

# Four tetrahedral vertices (unit vectors in local frame)
tetra_dirs = np.array([
    [ c,  s,  0],   # H1 direction (bonding)
    [ c, -s,  0],   # H2 direction (bonding)
    [-c,  0,  s],   # Lone pair 1 (back, +z)
    [-c,  0, -s],   # Lone pair 2 (back, -z)
])

print("Tetrahedral directions (local frame):")
print(f"  H1:  {tetra_dirs[0]}")
print(f"  H2:  {tetra_dirs[1]}")
print(f"  LP1: {tetra_dirs[2]}")
print(f"  LP2: {tetra_dirs[3]}")

# Verify angles
print("\nAngles between directions (should be 109.47° for tetrahedral):")
for i in range(4):
    for j in range(i+1, 4):
        angle = np.degrees(np.arccos(np.dot(tetra_dirs[i], tetra_dirs[j])))
        print(f"  {i}-{j}: {angle:.2f}°")

# sp3 orbital as linear combination with adjustable s:p ratio
# Standard sp3: c_s = 1/2, c_p = sqrt(3)/2  (s:p = 1 : sqrt(3))
# More p character -> more directional, density moves farther from nucleus
# c_s ranges from 0 (pure p) to 1 (pure s)

def hybrid_orbital(r_global, direction, c_s=0.5, zeta=2.0):
    """
    direction: unit vector [dx, dy, dz] pointing along the lobe
    c_s: s coefficient (0 to 1). c_p = sqrt(1 - c_s^2)
    Standard sp3: c_s = 0.5
    sp2: c_s = 1/sqrt(3) ≈ 0.577
    sp: c_s = 1/sqrt(2) ≈ 0.707
    """
    x, y, z = to_local(r_global)
    r = np.sqrt(x**2 + y**2 + z**2) + 1e-10
    s_orb = orbital_s(r, zeta)
    px = orbital_px(r, x, y, z, zeta)
    py = orbital_py(r, x, y, z, zeta)
    pz = orbital_pz(r, x, y, z, zeta)
    c_p = np.sqrt(1 - c_s**2)
    psi = c_s * s_orb + c_p * (direction[0]*px + direction[1]*py + direction[2]*pz)
    return psi

# Check actual H atom directions vs tetrahedral
H1_local = to_local(H1)
H2_local = to_local(H2)
print(f"\nActual H directions (normalized):")
print(f"  H1: {H1_local/np.linalg.norm(H1_local)}")
print(f"  H2: {H2_local/np.linalg.norm(H2_local)}")

# Lone pair directions are tetra_dirs[2] and tetra_dirs[3]
lone_pair_dirs = [tetra_dirs[2], tetra_dirs[3]]  # LP1 and LP2
bond_dirs = [tetra_dirs[0], tetra_dirs[1]]     # H1 and H2 bonding

# Grid parameters
L = 2.5
N = 200

# XY plane cross-section (z=0)
x_range = np.linspace(-L, L, N)
y_range = np.linspace(-L, L, N)
X_xy, Y_xy = np.meshgrid(x_range, y_range)

# YZ plane cross-section (x=0)
y_range_yz = np.linspace(-L, L, N)
z_range_yz = np.linspace(-L, L, N)
Y_yz, Z_yz = np.meshgrid(y_range_yz, z_range_yz)

# Vectorized compute density or amplitude
def compute_density_vec(plane, X, Y, Z, directions_list, c_s=0.5, zeta=2.2, clamp=False, power=1.0, return_amplitude=False):
    """Vectorized computation of |ψ|² or ψ for hybrid orbitals on a plane"""
    if plane == 'xy':
        # X,Y are 2D meshes; Z = 0
        r_global = np.stack([X, Y, np.zeros_like(X)], axis=-1)  # shape (N,N,3)
    elif plane == 'yz':
        r_global = np.stack([np.zeros_like(Y), Y, Z], axis=-1)
    elif plane == 'xz':
        r_global = np.stack([X, np.zeros_like(X), Z], axis=-1)

    # Convert to local frame: shape (N,N,3)
    dr = r_global - O  # broadcasting
    x_loc = np.einsum('ij,j->i', dr.reshape(-1, 3), x_axis).reshape(dr.shape[:2])
    y_loc = np.einsum('ij,j->i', dr.reshape(-1, 3), y_axis).reshape(dr.shape[:2])
    z_loc = np.einsum('ij,j->i', dr.reshape(-1, 3), z_axis).reshape(dr.shape[:2])

    r = np.sqrt(x_loc**2 + y_loc**2 + z_loc**2) + 1e-10
    s = orbital_s(r, zeta)

    c_p = np.sqrt(1 - c_s**2)
    result = np.zeros_like(r)
    for d in directions_list:
        # Angular part: dot product of direction with position
        angular = d[0]*x_loc + d[1]*y_loc + d[2]*z_loc
        if power == 1.0:
            # Standard: psi = c_s*s + c_p*angular*exp(-zeta*r)
            psi = c_s * s + c_p * angular * np.exp(-zeta * r)
        else:
            # Angular sharpening: psi = c_s*s + c_p*(angular/r)^power * r^power * exp(-zeta*r)
            # This gives (cos(theta))^power behavior
            cos_theta = angular / r
            psi = c_s * s + c_p * (cos_theta**power) * (r**power) * np.exp(-zeta * r)
        if clamp:
            psi = np.maximum(psi, 0)
        if return_amplitude:
            result += psi
        else:
            result += psi**2
    return result

# Compute densities using CLI parameters
c_p = np.sqrt(1 - c_s**2)
print(f"\nUsing: c_s={c_s:.3f}, c_p={c_p:.3f}, ratio s:p = 1:{c_p/c_s:.2f}")

dens_xy_lp1 = compute_density_vec('xy', X_xy, Y_xy, None, [lone_pair_dirs[0]], c_s=c_s, zeta=zeta, clamp=clamp, power=power, return_amplitude=amplitude)
dens_xy_total = compute_density_vec('xy', X_xy, Y_xy, None, lone_pair_dirs, c_s=c_s, zeta=zeta, clamp=clamp, power=power, return_amplitude=amplitude)

dens_yz_lp1 = compute_density_vec('yz', None, Y_yz, Z_yz, [lone_pair_dirs[0]], c_s=c_s, zeta=zeta, clamp=clamp, power=power, return_amplitude=amplitude)
dens_yz_total = compute_density_vec('yz', None, Y_yz, Z_yz, lone_pair_dirs, c_s=c_s, zeta=zeta, clamp=clamp, power=power, return_amplitude=amplitude)

# Helper for consistent contour levels
def plot_plane(ax, X, Y, data, title, xlabel, ylabel, atoms_x, atoms_y, is_amplitude=False):
    if is_amplitude:
        # Use diverging colormap for signed amplitude
        v_max = np.max(np.abs(data))
        im = ax.contourf(X, Y, data, levels=50, cmap='RdBu_r', vmin=-v_max, vmax=v_max)
        ax.contour(X, Y, data, levels=10, colors='black', alpha=0.3, linewidths=0.5)
        plt.colorbar(im, ax=ax, label='Amplitude ψ')
    else:
        # Use sequential colormap for density
        im = ax.contourf(X, Y, data, levels=50, cmap='viridis')
        ax.contour(X, Y, data, levels=10, colors='white', alpha=0.3, linewidths=0.5)
        plt.colorbar(im, ax=ax, label='Density |ψ|²')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_aspect('equal')
    ax.plot(atoms_x[0], atoms_y[0], 'ro', markersize=8, label='O')
    ax.plot(atoms_x[1:], atoms_y[1:], 'wo', markersize=5, label='H')
    ax.legend()
    return im

# Figure 1: Lone Pair 1 (individual)
fig1, axes1 = plt.subplots(1, 2, figsize=(12, 5))
plot_plane(axes1[0], X_xy, Y_xy, dens_xy_lp1,
           f'Lone Pair 1 (c_s={c_s:.2f}, dir: [{lone_pair_dirs[0][0]:.2f}, {lone_pair_dirs[0][1]:.2f}, {lone_pair_dirs[0][2]:.2f}]) - XY Plane', 'x (Å)', 'y (Å)',
           [O[0], H1[0], H2[0]], [O[1], H1[1], H2[1]], is_amplitude=amplitude)
plot_plane(axes1[1], Y_yz, Z_yz, dens_yz_lp1,
           f'Lone Pair 1 (c_s={c_s:.2f}, dir: [{lone_pair_dirs[0][0]:.2f}, {lone_pair_dirs[0][1]:.2f}, {lone_pair_dirs[0][2]:.2f}]) - YZ Plane', 'y (Å)', 'z (Å)',
           [O[1], H1[1], H2[1]], [O[2], H1[2], H2[2]], is_amplitude=amplitude)
plt.tight_layout()
plt.savefig(f'/home/prokophapala/git/FitHBonds/pyBall/{output_prefix}_lonepair1.png', dpi=150)

# Figure 2: Total (both lone pairs, incoherent sum |ψ1|²+|ψ2|²)
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
plot_plane(axes2[0], X_xy, Y_xy, dens_xy_total,
           f'Total Lone Pair Density (c_s={c_s:.2f}) - XY Plane', 'x (Å)', 'y (Å)',
           [O[0], H1[0], H2[0]], [O[1], H1[1], H2[1]], is_amplitude=amplitude)
plot_plane(axes2[1], Y_yz, Z_yz, dens_yz_total,
           f'Total Lone Pair Density (c_s={c_s:.2f}) - YZ Plane', 'y (Å)', 'z (Å)',
           [O[1], H1[1], H2[1]], [O[2], H1[2], H2[2]], is_amplitude=amplitude)
plt.tight_layout()
plt.savefig(f'/home/prokophapala/git/FitHBonds/pyBall/{output_prefix}_lonepairs_total.png', dpi=150)

plt.show()

print(f"\nPlots saved:")
print(f"  {output_prefix}_lonepair1.png    - individual orbital 1")
print(f"  {output_prefix}_lonepairs_total.png - incoherent sum of both")
