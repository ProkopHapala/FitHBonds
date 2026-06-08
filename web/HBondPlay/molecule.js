/* molecule.js - Molecular chemistry utilities for HBondPlay */

/* ---------- Atom Default Parameters ---------- */
const REQH_DEFAULTS = {
    'O':  {R:1.50, E:0.200, Q:-0.40, H:1.00},
    'H':  {R:1.00, E:0.050, Q:+0.20, H:0.00},
    'C':  {R:1.70, E:0.100, Q: 0.00, H:0.00},
    'N':  {R:1.50, E:0.200, Q:-0.40, H:1.00},
    'F':  {R:1.40, E:0.200, Q:-0.40, H:1.00},
    'Lp': {R:1.00, E:0.100, Q: 0.00, H:0.50}
};

/* ---------- QEq Parameters from ElementTypes.dat ---------- */
const QEQ_PARAMS = {
    'H':  { Eaff: -4.528,  Ehard: 13.8904 },
    'Li': { Eaff: -3.006,  Ehard: 4.772   },
    'C':  { Eaff: -5.343,  Ehard: 10.126  },
    'N':  { Eaff: -6.899,  Ehard: 11.760  },
    'O':  { Eaff: -8.741,  Ehard: 13.364  },
    'F':  { Eaff: -10.874, Ehard: 14.948  },
    'Na': { Eaff: -2.843,  Ehard: 4.592   },
    'Si': { Eaff: -4.168,  Ehard: 6.974   },
    'P':  { Eaff: -5.463,  Ehard: 8.000   },
    'S':  { Eaff: -6.928,  Ehard: 8.972   },
    'Se': { Eaff: -6.928,  Ehard: 8.972   },
    'Cl': { Eaff: -8.564,  Ehard: 9.892   },
    'Br': { Eaff: -8.564,  Ehard: 9.892   },
    'Ca': { Eaff: -2.843,  Ehard: 4.592   }
};

/* ---------- 3D Vector Math ---------- */
function normalizeVec(v) {
    const n = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if (n < 1e-8) return [0,0,1];
    return [v[0]/n, v[1]/n, v[2]/n];
}

function crossVec(a, b) {
    return [
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    ];
}

function scaleVec(v, s) {
    return [v[0]*s, v[1]*s, v[2]*s];
}

function addVec(a, b) {
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]];
}

function subVec(a, b) {
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
}

/* ---------- Electron Pair / Lone Pair Geometry ---------- */
function makeConfGeomJS(nb, npi, hs, Rlp) {
    // hs: array of bond direction vectors from central atom
    if (nb === 3) {
        const n = normalizeVec(crossVec(
            [hs[1][0]-hs[0][0], hs[1][1]-hs[0][1], hs[1][2]-hs[0][2]],
            [hs[2][0]-hs[0][0], hs[2][1]-hs[0][1], hs[2][2]-hs[0][2]]
        ));
        const dot = n[0]*(hs[0][0]+hs[1][0]+hs[2][0]) + n[1]*(hs[0][1]+hs[1][1]+hs[2][1]) + n[2]*(hs[0][2]+hs[1][2]+hs[2][2]);
        if (dot > 0) return [scaleVec(n, -Rlp)];
        return [scaleVec(n, Rlp)];
    } else if (nb === 2) {
        // Match C++ fromCrossSafe(hs[0], hs[1]):
        // b = normalize(cross(hs[0], hs[1]))
        // a = normalize(hs[1] - hs[0])
        // c = cross(b, a)
        const b = normalizeVec(crossVec(hs[0], hs[1]));
        const d = subVec(hs[1], hs[0]);
        const nd = Math.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        const a = nd > 1e-8 ? [d[0]/nd, d[1]/nd, d[2]/nd] : [1,0,0];
        const c = crossVec(b, a);
        if (npi === 0) {
            const cb = 0.81649658092;
            const cc = 0.57735026919;
            return [
                scaleVec(addVec(scaleVec(c, cc), scaleVec(b, cb)), Rlp),
                scaleVec(addVec(scaleVec(c, cc), scaleVec(b, -cb)), Rlp)
            ];
        } else {
            return [scaleVec(c, Rlp), scaleVec(b, Rlp)];
        }
    } else if (nb === 1) {
        const c = normalizeVec(hs[0]);
        let a = [1,0,0];
        if (Math.abs(c[0]) > 0.9) a = [0,1,0];
        let b = normalizeVec(crossVec(c, a));
        a = crossVec(b, c);
        if (npi === 0) {
            const ca = 0.81649658092;
            const cb = 0.47140452079;
            const cc = -0.33333333333;
            return [
                scaleVec(addVec(scaleVec(c, cc), scaleVec(b, cb*2)), Rlp),
                scaleVec(addVec(addVec(scaleVec(c, cc), scaleVec(b, -cb)), scaleVec(a, ca)), Rlp),
                scaleVec(addVec(addVec(scaleVec(c, cc), scaleVec(b, -cb)), scaleVec(a, -ca)), Rlp)
            ];
        } else if (npi === 1) {
            const ca = 0.86602540378;
            const cc = -0.5;
            return [
                scaleVec(addVec(scaleVec(c, cc), scaleVec(a, ca)), Rlp),
                scaleVec(addVec(scaleVec(c, cc), scaleVec(a, -ca)), Rlp),
                scaleVec(b, Rlp)
            ];
        } else {
            return [scaleVec(c, -Rlp), scaleVec(b, Rlp), scaleVec(a, Rlp)];
        }
    } else if (nb === 0) {
        const c = normalizeVec(hs[0]);
        let a = [1,0,0];
        if (Math.abs(c[0]) > 0.9) a = [0,1,0];
        let b = normalizeVec(crossVec(c, a));
        a = crossVec(b, c);
        if (npi === 0) {
            const ca = 0.81649658092;
            const cb = 0.47140452079;
            const cc = -0.33333333333;
            return [
                scaleVec(addVec(scaleVec(c, cc), scaleVec(b, cb*2)), Rlp),
                scaleVec(addVec(addVec(scaleVec(c, cc), scaleVec(b, -cb)), scaleVec(a, ca)), Rlp),
                scaleVec(addVec(addVec(scaleVec(c, cc), scaleVec(b, -cb)), scaleVec(a, -ca)), Rlp),
                scaleVec(c, Rlp)
            ];
        }
    }
    return [];
}

function generateEpairs(atoms, Rlp) {
    // Build neighbor list by distance
    const n = atoms.length;
    const neighs = Array(n).fill(null).map(() => []);
    const covalentR = {H:0.37, C:0.77, N:0.74, O:0.73, F:0.64};

    for (let i = 0; i < n; i++) {
        for (let j = i+1; j < n; j++) {
            const dx = atoms[i].x - atoms[j].x;
            const dy = atoms[i].y - atoms[j].y;
            const dz = atoms[i].z - atoms[j].z;
            const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
            const ri = covalentR[atoms[i].el] || 0.7;
            const rj = covalentR[atoms[j].el] || 0.7;
            if (dist < 1.3 * (ri + rj)) {
                neighs[i].push(j);
                neighs[j].push(i);
            }
        }
    }

    const result = atoms.map(a => ({...a}));
    const DEFAULT_EPAIRS = {C: 0, N: 1, O: 2, F: 3, Cl: 3, Br: 3, I: 3, H: 0};
    const HALOGENS = new Set(['F', 'Cl', 'Br', 'I']);

    for (let i = 0; i < n; i++) {
        const el = atoms[i].el;
        if (el === 'H') continue;
        const nb = neighs[i].length;
        let ne = DEFAULT_EPAIRS[el] || 0;

        // Special rule for halogens: 1 Lp (sigma hole) opposite to bond
        if (HALOGENS.has(el)) {
            ne = 1;
            if (nb >= 1 && neighs[i].length > 0) {
                const pos = [atoms[i].x, atoms[i].y, atoms[i].z];
                const j = neighs[i][0];
                const pj = [atoms[j].x, atoms[j].y, atoms[j].z];
                const bondDir = normalizeVec([pj[0]-pos[0], pj[1]-pos[1], pj[2]-pos[2]]);
                const lpDir = [-bondDir[0]*Rlp, -bondDir[1]*Rlp, -bondDir[2]*Rlp];
                result.push({
                    el: 'Lp', x: pos[0]+lpDir[0], y: pos[1]+lpDir[1], z: pos[2]+lpDir[2],
                    R: REQH_DEFAULTS.Lp.R, E: REQH_DEFAULTS.Lp.E,
                    Q: REQH_DEFAULTS.Lp.Q, H: REQH_DEFAULTS.Lp.H
                });
            }
            continue;
        }

        // Compute pi orbitals from octet rule: 4 = nb + npi + ne
        let npi = 4 - nb - ne;
        if (npi < 0) npi = 0;

        if (ne <= 0) continue;

        const pos = [atoms[i].x, atoms[i].y, atoms[i].z];
        const hs = neighs[i].map(j => {
            const pj = [atoms[j].x, atoms[j].y, atoms[j].z];
            return normalizeVec([pj[0]-pos[0], pj[1]-pos[1], pj[2]-pos[2]]);
        });

        if (hs.length === 0) hs.push([1,0,0]);

        const lpDirs = makeConfGeomJS(nb, npi, hs, Rlp);

        // Only take ne directions (makeConfGeom fills all 4 slots, we may need fewer)
        for (let k = 0; k < Math.min(ne, lpDirs.length); k++) {
            const d = lpDirs[k];
            result.push({
                el: 'Lp', x: pos[0]+d[0], y: pos[1]+d[1], z: pos[2]+d[2],
                R: REQH_DEFAULTS.Lp.R, E: REQH_DEFAULTS.Lp.E,
                Q: REQH_DEFAULTS.Lp.Q, H: REQH_DEFAULTS.Lp.H
            });
        }
    }

    return result;
}

/* ---------- XYZ Parser ---------- */
function parseXYZ(text) {
    const lines = text.trim().split(/\r?\n/);
    let idx = 0;
    const n = parseInt(lines[idx].trim());
    idx++;
    if (isNaN(n)) return null;
    // skip comment / lattice line if not atom-like
    while (idx < lines.length && !lines[idx].trim().match(/^[A-Za-z]/)) idx++;
    const atoms = [];
    for (let i = 0; i < n && idx < lines.length; i++, idx++) {
        const parts = lines[idx].trim().split(/\s+/);
        if (parts.length < 4) { i--; continue; }
        atoms.push({
            el: parts[0],
            x: parseFloat(parts[1]),
            y: parseFloat(parts[2]),
            z: parseFloat(parts[3])
        });
    }
    return atoms;
}

/* ---------- QEq Charge Equilibration (ported from cpp/molecular/QEq.h) ---------- */
function makeCoulombMatrix(n, ps) {
    // ps: array of {x,y,z} objects or arrays
    const J = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
        const pi = ps[i];
        J[i * n + i] = 0.0;
        for (let j = i + 1; j < n; j++) {
            const pj = ps[j];
            const dx = pj.x !== undefined ? pj.x - pi.x : pj[0] - pi[0];
            const dy = pj.y !== undefined ? pj.y - pi.y : pj[1] - pi[1];
            const dz = pj.z !== undefined ? pj.z - pi.z : pj[2] - pi[2];
            const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
            const Jij = 14.3996448915 / (1.0 + r);
            J[i * n + j] = Jij;
            J[j * n + i] = Jij;
        }
    }
    return J;
}

class QEq {
    constructor() {
        this.n = 0;
        this.J = null;
        this.qs = null;
        this.fqs = null;
        this.vqs = null;
        this.affins = null;
        this.hards = null;
        this.constrain = null;
        this.Qtarget = 0.0;
        this.Qtot = 0.0;
    }

    realloc(n_) {
        this.n = n_;
        this.qs = new Float64Array(n_);
        this.fqs = new Float64Array(n_);
        this.vqs = new Float64Array(n_);
        this.hards = new Float64Array(n_);
        this.affins = new Float64Array(n_);
        this.constrain = new Array(n_).fill(false);
    }

    init() {
        for (let i = 0; i < this.n; i++) {
            this.qs[i] = 0.0;
            this.fqs[i] = 0.0;
            this.vqs[i] = 0.0;
        }
    }

    constrainTypes(atypes, iconst) {
        for (let i = 0; i < this.n; i++) {
            if (atypes[i] === iconst) this.constrain[i] = true;
        }
    }

    getQvars() {
        let err2 = 0.0;
        let fqtot = 0.0;
        let nsum = 0;
        for (let i = 0; i < this.n; i++) {
            if (this.constrain[i]) continue;
            const qi = this.qs[i];
            let fq = this.affins[i] + this.hards[i] * qi;
            for (let j = 0; j < this.n; j++) {
                fq += this.J[i * this.n + j] * this.qs[j];
            }
            fqtot += fq;
            this.fqs[i] = fq;
            nsum++;
        }
        const dfqtot = fqtot / nsum;
        for (let i = 0; i < this.n; i++) {
            if (this.constrain[i]) continue;
            this.fqs[i] -= dfqtot;
            err2 += this.fqs[i] * this.fqs[i];
        }
        return err2;
    }

    moveMDdamp(dt, damping) {
        this.Qtot = 0.0;
        const damp = 1.0 - damping;
        let nsum = 0;
        let cvf = 0.0;
        for (let i = 0; i < this.n; i++) {
            let qi = this.qs[i];
            if (this.constrain[i]) {
                this.Qtot += qi;
                continue;
            }
            const fqi = this.fqs[i];
            let vqi = this.vqs[i];
            vqi *= damp;
            vqi -= fqi * dt;
            qi += vqi * dt;
            this.vqs[i] = vqi;
            this.qs[i] = qi;
            this.Qtot += qi;
            cvf += fqi * vqi;
            nsum++;
        }
        const dQ = (this.Qtarget - this.Qtot) / nsum;
        for (let i = 0; i < this.n; i++) {
            if (this.constrain[i]) continue;
            this.qs[i] += dQ;
        }
        return this.Qtot;
    }

    relaxChargeMD(ps, nsteps = 1000, Fconf = 1e-2, dt = 0.1, damp = 0.0, bVerbose = false, bInit = true) {
        const F2conf = Fconf * Fconf;
        this.J = makeCoulombMatrix(this.n, ps);
        let F2 = 1.0;
        if (bInit) {
            this.init();
        } else {
            for (let i = 0; i < this.n; i++) this.qs[i] *= -1.0;
        }
        for (let itr = 0; itr < nsteps; itr++) {
            F2 = this.getQvars();
            if (F2 < F2conf) break;
            if (bVerbose) console.log(`QEq.relaxChargeMD()[${itr}] |F|=${Math.sqrt(F2)}`);
            const cvf = this.moveMDdamp(dt, damp);
            if (cvf < 0) {
                if (bVerbose) console.log(`QEq.relaxChargeMD()[${itr}](cvf(${cvf})<0) => v[:]=0`);
                for (let i = 0; i < this.n; i++) this.vqs[i] = 0.0;
            }
        }
        let Qtot = 0.0;
        for (let i = 0; i < this.n; i++) {
            this.qs[i] *= -1.0;
            Qtot += this.qs[i];
        }
        if (Math.abs(Qtot - this.Qtarget) > 1e-6) {
            console.warn(`QEq.relaxChargeMD() WARNING: Qtot(${Qtot})-Qtarget(${this.Qtarget})=${Qtot - this.Qtarget}`);
        }
        this.J = null;
        return F2;
    }
}

function runQEqOnFragment(atoms, Qtarget = 0.0) {
    // Filter out lone pairs (Lp) from QEq - they are not real atoms
    const atomIndices = [];
    const positions = [];
    for (let i = 0; i < atoms.length; i++) {
        if (atoms[i].el === 'Lp') continue;
        atomIndices.push(i);
        positions.push(atoms[i]);
    }
    const n = atomIndices.length;
    if (n === 0) return;

    const qeq = new QEq();
    qeq.realloc(n);
    qeq.Qtarget = Qtarget;

    for (let i = 0; i < n; i++) {
        const el = atoms[atomIndices[i]].el;
        const p = QEQ_PARAMS[el];
        if (p) {
            qeq.affins[i] = p.Eaff;
            qeq.hards[i] = p.Ehard;
        } else {
            console.warn(`No QEq params for element ${el}, using defaults`);
            qeq.affins[i] = -5.0;
            qeq.hards[i] = 10.0;
        }
    }

    qeq.relaxChargeMD(positions, 2000, 1e-4, 0.1, 0.05, false, true);

    // Write back charges (proton convention)
    for (let i = 0; i < n; i++) {
        atoms[atomIndices[i]].Q = qeq.qs[i];
    }
}

/* ---------- Molecule Class ---------- */
class Molecule {
    constructor(name, atoms, type) {
        this.name = name;
        this.atoms = atoms; // {el, x, y, z, R, E, Q, H}
        this.type = type;
        this.n = atoms.length;
    }

    center() {
        const ox = this.atoms[0].x, oy = this.atoms[0].y, oz = this.atoms[0].z;
        this.atoms.forEach(a => { a.x -= ox; a.y -= oy; a.z -= oz; });
    }

    attachDefaults() {
        this.atoms.forEach(a => {
            const d = REQH_DEFAULTS[a.el] || {R:1.0, E:0.1, Q:0.0, H:0.0};
            a.R = d.R; a.E = d.E; a.Q = d.Q; a.H = d.H;
        });
    }

    addEpairs(Rlp) {
        this.atoms = generateEpairs(this.atoms, Rlp);
        this.n = this.atoms.length;
        let lpTotalQ = 0;
        for (let i = 0; i < this.n; i++) {
            if (this.atoms[i].el === 'Lp') lpTotalQ += this.atoms[i].Q;
        }
        if (lpTotalQ !== 0) {
            const nHosts = this.atoms.filter(a => a.el !== 'Lp').length;
            if (nHosts > 0) {
                const adjustment = -lpTotalQ / nHosts;
                this.atoms.forEach(a => { if (a.el !== 'Lp') a.Q += adjustment; });
            }
        }
    }

    rotateToX() {
        let far = this.atoms[1] || this.atoms[0];
        let farD = far.x*far.x + far.y*far.y + far.z*far.z;
        for (let i = 2; i < this.n; i++) {
            const d = this.atoms[i].x*this.atoms[i].x + this.atoms[i].y*this.atoms[i].y + this.atoms[i].z*this.atoms[i].z;
            if (d > farD) { farD = d; far = this.atoms[i]; }
        }
        const v = [far.x, far.y, far.z];
        const vl = Math.sqrt(farD);
        if (vl > 0.001) {
            v[0] /= vl; v[1] /= vl; v[2] /= vl;
        } else {
            v[0] = 1; v[1] = 0; v[2] = 0;
        }
        const axis = [0, v[2], -v[1]];
        const al = Math.sqrt(axis[1]*axis[1] + axis[2]*axis[2]);
        let angle = Math.acos(v[0]);
        if (al < 0.001) {
            if (v[0] < 0) {
                this.atoms.forEach(a => { a.x = -a.x; a.y = -a.y; });
            }
        } else {
            axis[0] /= al; axis[1] /= al; axis[2] /= al;
            const c = Math.cos(angle), s = Math.sin(angle);
            const ux=axis[0], uy=axis[1], uz=axis[2];
            const rot = (p) => {
                const dot = ux*p.x + uy*p.y + uz*p.z;
                const cx = uy*p.z - uz*p.y;
                const cy = uz*p.x - ux*p.z;
                const cz = ux*p.y - uy*p.x;
                return {
                    ...p,
                    x: p.x*c + cx*s + ux*dot*(1-c),
                    y: p.y*c + cy*s + uy*dot*(1-c),
                    z: p.z*c + cz*s + uz*dot*(1-c)
                };
            };
            this.atoms = this.atoms.map(rot);
        }
    }

    // Evaluate interaction energy between two molecules at a given placement.
    // mol1 is fixed, mol2 is placed with its local origin at 'center',
    // with its local x-axis along 'dir', y along 'perp', z along 'normal'.
    static energyAt(mol1, mol2, center, dir, perp, normal, b, mask) {
        let E = 0;
        for (let i = 0; i < mol1.n; i++) {
            const ai = mol1.atoms[i];
            for (let j = 0; j < mol2.n; j++) {
                const aj = mol2.atoms[j];
                const px = center[0] + dir[0]*aj.x + perp[0]*aj.y + normal[0]*aj.z;
                const py = center[1] + dir[1]*aj.x + perp[1]*aj.y + normal[1]*aj.z;
                const pz = center[2] + dir[2]*aj.x + perp[2]*aj.y + normal[2]*aj.z;
                const dx = ai.x - px, dy = ai.y - py, dz = ai.z - pz;
                let r = Math.sqrt(dx*dx + dy*dy + dz*dz);
                if (r < 0.05) r = 0.05;
                const Rij = ai.R + aj.R;
                const Eij = Math.sqrt(Math.max(0, ai.E * aj.E));
                const dr = r - Rij;
                const eb = Math.exp(-b * dr);
                const eb2 = eb*eb;
                if (mask[0]) E += Eij * eb2;
                if (mask[1]) E -= 2.0 * Eij * eb;
                if (mask[2]) E += 332.0637 * ai.Q * aj.Q / r;
            }
        }
        return E;
    }

    // 1D line profile along the plot y-axis through the pivot.
    // In polar mode:  center angle = (aMin+aMax)/2, distance varies.
    // In cartesian:    x = 0, y varies.
    static lineProfile(mol1, mol2, pivotIdx, mode, aMin, aMax, dMin, dMax, cartRange, b, mask, nSamples = 64) {
        const profile = [];
        for (let iy = 0; iy <= nSamples; iy++) {
            let dir, perp, normal = [0,0,1], cx, cy, cz = 0;
            if (mode === 'polar') {
                const ang = (aMin + aMax) * 0.5;
                const dist = dMin + (dMax - dMin) * iy / nSamples;
                const ca = Math.cos(ang), sa = Math.sin(ang);
                dir = [sa, ca, 0];
                perp = [-ca, sa, 0];
                cx = dist * dir[0];
                cy = dist * dir[1];
            } else {
                const worldX = 0;
                const worldY = -cartRange + (2 * cartRange) * iy / nSamples;
                cx = worldX; cy = worldY;
                const toPixel = [worldX, worldY, 0];
                const tpLen = Math.sqrt(toPixel[0]*toPixel[0] + toPixel[1]*toPixel[1]);
                if (tpLen < 0.001) {
                    dir = [1, 0, 0];
                } else {
                    dir = [toPixel[0]/tpLen, toPixel[1]/tpLen, 0];
                }
                perp = [-dir[1], dir[0], 0];
            }
            const E = Molecule.energyAt(mol1, mol2, [cx, cy, cz], dir, perp, normal, b, mask);
            const yVal = mode === 'polar'
                ? (dMin + (dMax - dMin) * iy / nSamples)
                : (-cartRange + (2 * cartRange) * iy / nSamples);
            profile.push({ y: yVal, E });
        }
        return profile;
    }
}
