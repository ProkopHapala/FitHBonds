"""
Polynomial Basis for Atomic Orbitals
Based on doc/ideas/Polynominal_AOs.md

This script uses sympy to compute exact Gram-Schmidt orthogonalization
coefficients with the proper r^2 Jacobian weighting for radial functions.

The inner product for radial functions in 3D includes the spherical
volume element: <f|g> = integral_0^Rc f(r) g(r) r^2 dr

In dimensionless variable x = r/Rc, this becomes:
    <P_i|P_j> = Rc^3 * integral_0^1 (1-x)^(i+j) * x^2 dx

The integral evaluates to the Beta function:
    I(k) = 2 / ((k+1)(k+2)(k+3))   where k = i+j
"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, integrate, sqrt, Rational, pprint, simplify

# ============================================================
# Sympy exact computation
# ============================================================

x_sym = symbols('x', real=True, positive=True)

def inner_product_exact(i, j):
    """Exact symbolic inner product <P_i|P_j> with r^2 weighting."""
    # P_i(x) = (1-x)^i,  P_j(x) = (1-x)^j
    # Inner product = integral_0^1 (1-x)^(i+j) * x^2 dx
    integrand = (1 - x_sym)**(i + j) * x_sym**2
    val = integrate(integrand, (x_sym, 0, 1))
    return val


def gram_schmidt_exact(powers, top_down=True):
    """
    Exact symbolic Gram-Schmidt orthogonalization.
    
    Returns exact rational coefficients and symbolic norms.
    """
    if top_down:
        order = sorted(powers, reverse=True)
    else:
        order = sorted(powers)
    
    n = len(order)
    
    # Precompute all inner products
    S = {}
    for p1 in powers:
        for p2 in powers:
            S[(p1, p2)] = inner_product_exact(p1, p2)
    
    coeffs_exact = []
    
    for m in range(n):
        # Start with seed polynomial
        c = {p: 0 for p in powers}
        c[order[m]] = Rational(1)
        
        # Orthogonalize against all previous functions
        for k in range(m):
            # Previous function: psi_k = sum_p coeffs_exact[k][p] * P_p
            # Compute <P_order[m] | psi_k>
            num = Rational(0)
            for p in powers:
                if coeffs_exact[k][p] != 0:
                    num += coeffs_exact[k][p] * S[(order[m], p)]
            
            # Compute <psi_k | psi_k>
            den = Rational(0)
            for p1 in powers:
                for p2 in powers:
                    if coeffs_exact[k][p1] != 0 and coeffs_exact[k][p2] != 0:
                        den += coeffs_exact[k][p1] * coeffs_exact[k][p2] * S[(p1, p2)]
            
            alpha = num / den
            for p in powers:
                c[p] -= alpha * coeffs_exact[k][p]
        
        coeffs_exact.append(c)
    
    # Normalize each function
    coeffs_norm = []
    for m in range(n):
        norm_sq = Rational(0)
        for p1 in powers:
            for p2 in powers:
                if coeffs_exact[m][p1] != 0 and coeffs_exact[m][p2] != 0:
                    norm_sq += coeffs_exact[m][p1] * coeffs_exact[m][p2] * S[(p1, p2)]
        
        norm = sqrt(norm_sq)
        coeffs_norm.append({p: simplify(coeffs_exact[m][p] / norm) for p in powers})
    
    return coeffs_norm, order


# ============================================================
# Numerical evaluation for plotting
# ============================================================

def P(n, x):
    """Raw polynomial P_n(x) = (1-x)^n for x in [0,1], zero beyond."""
    return np.where(x <= 1.0, (1.0 - x) ** n, 0.0)


def evaluate_combination(coeffs_dict, powers, x):
    """Evaluate a linear combination of raw polynomials from sympy coefficients."""
    result = np.zeros_like(x, dtype=float)
    for p in powers:
        c = float(coeffs_dict[p])
        if abs(c) > 1e-15:
            result += c * P(p, x)
    return result


# ============================================================
# Plotting
# ============================================================

def plot_sequence(powers, title_suffix, ax_raw, ax_ortho_td, ax_ortho_bu):
    """Plot raw and orthogonalized polynomials for one sequence."""
    x = np.linspace(0, 1.2, 500)
    
    # --- Raw powers ---
    for p in powers:
        y = P(p, x)
        ax_raw.plot(x, y, label=f"$P_{{{p}}}$")
    ax_raw.set_title(f"Raw: {title_suffix}")
    ax_raw.set_xlabel("$x = r/R_c$")
    ax_raw.set_ylabel("$P_n(x)$")
    ax_raw.legend(loc='upper right')
    ax_raw.axvline(1.0, color='k', linestyle='--', alpha=0.3)
    ax_raw.set_xlim(0, 1.2)
    ax_raw.grid(True, alpha=0.3)
    
    # --- Top-down orthogonalization ---
    coeffs_td, order_td = gram_schmidt_exact(powers, top_down=True)
    colors = plt.cm.tab10(np.linspace(0, 0.9, len(powers)))
    
    for i, (c, p_seed) in enumerate(zip(coeffs_td, order_td)):
        y = evaluate_combination(c, powers, x)
        label = f"$\\psi_{{{i}}}$ seed $P_{{{p_seed}}}$"
        ax_ortho_td.plot(x, y, label=label, color=colors[i], linewidth=1.5)
    
    ax_ortho_td.set_title(f"Top-down: {title_suffix}")
    ax_ortho_td.set_xlabel("$x = r/R_c$")
    ax_ortho_td.set_ylabel("$\\psi_i(x)$")
    ax_ortho_td.legend(loc='upper right', fontsize='small')
    ax_ortho_td.axvline(1.0, color='k', linestyle='--', alpha=0.3)
    ax_ortho_td.set_xlim(0, 1.2)
    ax_ortho_td.grid(True, alpha=0.3)
    
    # Print exact coefficients
    print(f"\n{'='*60}")
    print(f"Top-down: {title_suffix}   (seed order: {order_td})")
    print(f"{'='*60}")
    for i, (c, p_seed) in enumerate(zip(coeffs_td, order_td)):
        terms = []
        for p in powers:
            coeff = c[p]
            if coeff != 0:
                terms.append(f"({coeff})*P{p}")
        print(f"psi_{i} (seed P_{p_seed}): {' + '.join(terms)}")
    
    # --- Bottom-up orthogonalization ---
    coeffs_bu, order_bu = gram_schmidt_exact(powers, top_down=False)
    
    for i, (c, p_seed) in enumerate(zip(coeffs_bu, order_bu)):
        y = evaluate_combination(c, powers, x)
        label = f"$\\psi_{{{i}}}$ seed $P_{{{p_seed}}}$"
        ax_ortho_bu.plot(x, y, label=label, color=colors[i], linewidth=1.5)
    
    ax_ortho_bu.set_title(f"Bottom-up: {title_suffix}")
    ax_ortho_bu.set_xlabel("$x = r/R_c$")
    ax_ortho_bu.set_ylabel("$\\psi_i(x)$")
    ax_ortho_bu.legend(loc='upper right', fontsize='small')
    ax_ortho_bu.axvline(1.0, color='k', linestyle='--', alpha=0.3)
    ax_ortho_bu.set_xlim(0, 1.2)
    ax_ortho_bu.grid(True, alpha=0.3)
    
    print(f"\n{'='*60}")
    print(f"Bottom-up: {title_suffix}   (seed order: {order_bu})")
    print(f"{'='*60}")
    for i, (c, p_seed) in enumerate(zip(coeffs_bu, order_bu)):
        terms = []
        for p in powers:
            coeff = c[p]
            if coeff != 0:
                terms.append(f"({coeff})*P{p}")
        print(f"psi_{i} (seed P_{p_seed}): {' + '.join(terms)}")


# ============================================================
# Main execution
# ============================================================

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Sequence 1: p2, p4, p6, p8
powers_2468 = [2, 4, 6, 8]
plot_sequence(powers_2468, "$P_2, P_4, P_6, P_8$",
                axes[0, 0], axes[0, 1], axes[0, 2])

# Sequence 2: p2, p4, p8, p16
powers_24816 = [2, 4, 8, 16]
plot_sequence(powers_24816, "$P_2, P_4, P_8, P_{{16}}$",
                axes[1, 0], axes[1, 1], axes[1, 2])

plt.tight_layout()
plt.savefig("polynomial_basis.png", dpi=150)
plt.savefig("polynomial_basis.svg")
print("\nSaved: polynomial_basis.png, polynomial_basis.svg")

# Also create separate focused plots for raw powers
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
x = np.linspace(0, 1.2, 500)

for p in [2, 4, 6, 8]:
    axes2[0].plot(x, P(p, x), label=f"$P_{{{p}}}$")
axes2[0].set_title("Raw powers: $P_2, P_4, P_6, P_8$")
axes2[0].set_xlabel("$x = r/R_c$")
axes2[0].set_ylabel("$P_n(x)$")
axes2[0].legend(loc='upper right')
axes2[0].axvline(1.0, color='k', linestyle='--', alpha=0.3)
axes2[0].set_xlim(0, 1.2)
axes2[0].grid(True, alpha=0.3)

for p in [2, 4, 8, 16]:
    axes2[1].plot(x, P(p, x), label=f"$P_{{{p}}}$")
axes2[1].set_title("Raw powers: $P_2, P_4, P_8, P_{{16}}$")
axes2[1].set_xlabel("$x = r/R_c$")
axes2[1].set_ylabel("$P_n(x)$")
axes2[1].legend(loc='upper right')
axes2[1].axvline(1.0, color='k', linestyle='--', alpha=0.3)
axes2[1].set_xlim(0, 1.2)
axes2[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("raw_powers.png", dpi=150)
plt.savefig("raw_powers.svg")
print("Saved: raw_powers.png, raw_powers.svg")
