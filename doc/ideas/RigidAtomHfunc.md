https://gemini.google.com/share/28851650a5ee

---

# USER

Rigid atom hydrogen bonds.

I would like to make hydrogen bond corrections by adding angular functions to atoms but instead of computing orientation of coordinate system on the fly from bonds to neighbors like in other forcefield each atom should carry it's own reference frame (quaternion) as a rigid body. In this reference frame it can have few angular functions maybe one or two referring direction of electron pairs or signal holes. But these are not charges decaying like 1\r. These are electron cloud polarization functions like Fukui function (represent willingness to give or take electron, electrophilic or nucleophilic attack susceptibility). They decaying exponential from atom. Imagine e.g. H2O with two electron pairs or hydrogen or bromine with sigma hole. These can be like SP3,SP2,SP1 orbitals made by combination of s and p function. But we represent only the density (square) of that function.

---

# Gemini

This is a very elegant and computationally efficient way to handle anisotropic molecular interactions. By decoupling the local orientation from the dynamic neighbor list and tracking it via an explicit quaternion, you completely bypass the nasty, non-continuous derivative jumps that happen in traditional forcefields when neighbor definitions switch or bonds stretch.

Furthermore, using exponential decay ($\sim e^{-\alpha r}$) instead of electrostatic $1/r$ is physically much closer to the true nature of hydrogen bonding, which relies heavily on **wavefunction overlap, charge transfer, and Pauli repulsion** at short range.

Here is a blueprint for how you can mathematically frame and implement this "Rigid Atom with Exponential Angular Density" model.

## 1. Local Frame and Angular Representation

Each atom $i$ carries a position $\vec{R}_i$ and an orientation quaternion $\hat{q}_i$. This quaternion defines a rotation matrix $\mathbf{R}(\hat{q}_i)$ that transforms vectors from the local atomic coordinate system to the global simulation box.

Instead of full orbital wavefunctions, you are tracking the *density* (probability distribution) of the interaction zones (lone pairs, $\sigma$-holes). We can represent these as a sum of angular functions fixed in the atom's local frame, multiplied by a radial exponential decay.

### The Angular Density Function

For a specific interaction site $\alpha$ on atom $i$, the local density contribution $F_{i,\alpha}$ in the direction of a unit vector $\vec{u}$ can be written using a combination of spherical harmonics (or simple dot products aligned with local hybrid-like directions $\vec{d}_{i,\alpha}$):

$$F_{i,\alpha}(\vec{u}) = c_{i,\alpha} \cdot \max(0, \vec{u} \cdot \vec{d}_{i,\alpha})^{n}$$

* **$\vec{d}_{i,\alpha}$**: A fixed unit vector in the atom's local frame pointing toward the lone pair or $\sigma$-hole (e.g., for $sp^3$ water, these would point toward the tetrahedral vertices).
* **$n$**: A power controlling the angular sharpness (directional lock) of the interaction.
* **$c_{i,\alpha}$**: A coefficient representing the character (e.g., positive for an electrophilic $\sigma$-hole / Lewis acid, negative for a nucleophilic lone pair / Lewis base).

## 2. Interaction Energy Formulation

When atom $i$ (donor/acceptor) interacts with atom $j$, the interaction depends on the distance $r_{ij} = |\vec{R}_j - \vec{R}_i|$ and the relative orientations.

Let $\vec{u}_{ij} = (\vec{R}_j - \vec{R}_i) / r_{ij}$ be the global unit vector pointing from $i$ to $j$. To evaluate the angular functions, we project this vector into the respective local frames of both atoms:

$$\vec{u}_{i}^{\text{local}} = \mathbf{R}(\hat{q}_i)^T \vec{u}_{ij}$$

$$\vec{u}_{j}^{\text{local}} = \mathbf{R}(\hat{q}_j)^T (-\vec{u}_{ij})$$

### The Pair Potential

The directional correction energy $E_{\text{HB}}$ between atom $i$ and atom $j$ can be modeled as the overlap susceptibility of their respective functions, decaying exponentially with distance:

$$E_{\text{HB}}(i,j) = \epsilon_{ij} \cdot \left( \sum_{\alpha \in i} F_{i,\alpha}(\vec{u}_{i}^{\text{local}}) \right) \cdot \left( \sum_{\beta \in j} F_{j,\beta}(\vec{u}_{j}^{\text{local}}) \right) \cdot e^{-\beta_{ij} (r_{ij} - r_0)}$$

* **$\epsilon_{ij}$**: Interaction strength scaling factor. If you match a positive Fukui-like function ($\sigma$-hole) with a negative one (lone pair), $\epsilon_{ij}$ would be negative to ensure attraction.
* **$\beta_{ij}$**: The exponential decay constant (typically chosen to match the overlap of atomic densities, e.g., around $2.0 \text{ to } 4.0 \, \mathring{\text{A}}^{-1}$).
* **$r_0$**: A reference equilibrium distance.

> **Note on Pauli Repulsion:** Because this function decays exponentially, you can seamlessly pair it with an exponential repulsion term ($\sim e^{-2\beta r}$) to form a Buckingham-like potential for the hydrogen bond, completely eliminating the need for Lennard-Jones $1/r^{12}$ steepness in the hydrogen bond channel.

## 3. Forces and Torques

Since the atoms behave as rigid bodies, the energy gradients must be mapped to both **translational forces** and **rotational torques**.

### Translational Force

The force on atom $i$ is the negative gradient with respect to positions:


$$\vec{F}_i = -\nabla_{\vec{R}_i} E_{\text{HB}}$$


Because $\vec{u}_{ij}$ depends on positions, $\nabla$ will hit both the radial exponential part and the angular unit vectors.

### Torques (The Rigid Body Update)

This is where the quaternion reference frame shines. The torque $\vec{\tau}_i$ acting on atom $i$ in the global frame is derived from how the energy changes when the local frame rotates:

$$\vec{\tau}_i = - \vec{u}_{ij} \times \nabla_{\vec{u}_{ij}} E_{\text{HB}}$$

In practice, for a single interaction pair, the torque generated by the angular function on atom $i$ is:

$$\vec{\tau}_i = \epsilon_{ij} \cdot e^{-\beta_{ij}(r_{ij}-r_0)} \cdot \left(\sum_{\beta \in j} F_j \right) \cdot \left[ \mathbf{R}(\hat{q}_i) \left( \vec{u}_{i}^{\text{local}} \times \nabla_{\vec{u}_{i}^{\text{local}}} F_i \right) \right]$$

You evaluate the derivative of the angular shape in the easy local frame, take the cross product with the local look-vector, and then rotate that torque vector back to the global frame using $\mathbf{R}(\hat{q}_i)$.

## 4. Mapping Specific Chemical Motifs

Because you are using local densities rather than explicit wavefunctions, you can construct standard hybrid profiles easily:

| Motif | Local Vector Alignment ($\vec{d}$) | Sign ($c$) | Physical Meaning |
| --- | --- | --- | --- |
| **Water Oxygen ($sp^3$)** | Two vectors pointing at $\approx 109.5^\circ$ | Negative | Lone pairs (Nucleophilic) |
| **Water Hydrogen** | One vector pointing along the O-H axis | Positive | Electron deficiency (Electrophilic) |
| **Bromine / Iodine ($\sigma$-hole)** | One vector pointing directly along the C-X bond axis | Positive | Electropositive crown on the tip of the halogen |
| **Carbonyl Oxygen ($sp^2$)** | Two vectors in the molecular plane at $\approx 120^\circ$ to the C=O bond | Negative | $sp^2$ lone pairs |

## Implementation Advantages

1. **No Discontinuities:** Traditional forcefields redefine local axes based on $3$-body or $4$-body connectivity. If a hydrogen bond geometry distorts severely, or if you are dealing with a molten/reactive system where bonds break, those definitions fail. Quaternions don't care about connectivity; they evolve smoothly via rigid-body equations of motion.
2. **GPU-Friendly:** Evaluating a quaternion rotation matrix and a few dot products requires minimal registers and zero branch hopping, making this highly parallelizable for OpenCL/CUDA or high-performance Rust implementations.

How do you plan to handle the time propagation of the quaternions? Are you leaning toward a standard integration scheme like NoSQUISH or a direct angular velocity update?