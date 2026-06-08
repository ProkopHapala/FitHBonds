/* HBondPlay app.js */

const MOLECULES = {
    'H2O':  { file: 'data/xyz/H2O.xyz',  type: 'frag1', addLP: false },
    'CH2O': { file: 'data/xyz/CH2O.xyz', type: 'frag1', addLP: false },
    'CH2NH':{ file: 'data/xyz/CH2NH.xyz',type: 'frag1', addLP: false },
    'HCN':  { file: 'data/xyz/HCN.xyz',  type: 'frag2', addLP: false },
    'HF':   { file: 'data/xyz/HF.xyz',   type: 'frag2', addLP: false }
};

const COLORMAPS = {
    'viridis':  null,
    'magma':    null,
    'inferno':  null,
    'plasma':   null,
    'coolwarm': null,
    'seismic':  null
};

let gl, program, locs;
let current = {
    frag1: null, frag2: null,
    pivot1: 0,
    debugClick: null
};

async function loadMolecule(name) {
    const info = MOLECULES[name];
    const resp = await fetch(info.file);
    const text = await resp.text();
    let atoms = parseXYZ(text);
    if (!atoms) return null;

    const mol = new Molecule(name, atoms, info.type);
    mol.center();
    mol.attachDefaults();

    if (info.type === 'frag1' && document.getElementById('addEpair')?.checked) {
        const Rlp = parseFloat(document.getElementById('lpDist').value);
        mol.addEpairs(Rlp);
    }
    if (info.type === 'frag2') {
        mol.rotateToX();
    }
    return mol;
}

/* ---------- Colormaps ---------- */
function buildColormaps() {
    COLORMAPS.viridis = genMap((t) => {
        const c = [[0.267,0.004,0.329],[0.283,0.141,0.458],[0.253,0.265,0.529],[0.163,0.471,0.558],[0.135,0.620,0.537],[0.271,0.749,0.407],[0.627,0.855,0.278],[0.993,0.906,0.144]];
        return spline8(t, c);
    });
    COLORMAPS.magma = genMap((t) => {
        const c = [[0.001,0.000,0.013],[0.183,0.070,0.338],[0.504,0.154,0.471],[0.752,0.252,0.471],[0.916,0.371,0.424],[0.978,0.518,0.392],[0.989,0.705,0.518],[0.988,0.937,0.565]];
        return spline8(t, c);
    });
    COLORMAPS.inferno = genMap((t) => {
        const c = [[0.001,0.000,0.014],[0.169,0.026,0.278],[0.482,0.074,0.361],[0.783,0.179,0.304],[0.986,0.350,0.196],[0.974,0.553,0.174],[0.940,0.752,0.287],[0.988,0.998,0.645]];
        return spline8(t, c);
    });
    COLORMAPS.plasma = genMap((t) => {
        const c = [[0.050,0.027,0.529],[0.321,0.060,0.644],[0.533,0.147,0.615],[0.723,0.251,0.535],[0.890,0.371,0.432],[0.997,0.518,0.302],[0.988,0.718,0.149],[0.940,0.975,0.131]];
        return spline8(t, c);
    });
    COLORMAPS.coolwarm = genMap((t) => {
        const c = [[0.230,0.299,0.754],[0.436,0.650,0.889],[0.865,0.865,0.865],[0.956,0.647,0.509],[0.706,0.016,0.150]];
        return spline5(t, c);
    });
    COLORMAPS.seismic = genMap((t) => {
        const c = [[0.000,0.000,0.300],[0.000,0.000,1.000],[1.000,1.000,1.000],[1.000,0.000,0.000],[0.500,0.000,0.000]];
        return spline5(t, c);
    });
}

function genMap(fn) {
    const arr = new Uint8Array(256*3);
    for (let i = 0; i < 256; i++) {
        const t = i / 255;
        const rgb = fn(t);
        arr[i*3]   = Math.max(0, Math.min(255, rgb[0]*255));
        arr[i*3+1] = Math.max(0, Math.min(255, rgb[1]*255));
        arr[i*3+2] = Math.max(0, Math.min(255, rgb[2]*255));
    }
    return arr;
}

function spline8(t, cp) {
    const s = t * 7;
    const i = Math.min(6, Math.floor(s));
    const u = s - i;
    const a = cp[i], b = cp[i+1];
    return [a[0]+(b[0]-a[0])*u, a[1]+(b[1]-a[1])*u, a[2]+(b[2]-a[2])*u];
}
function spline5(t, cp) {
    const s = t * 4;
    const i = Math.min(3, Math.floor(s));
    const u = s - i;
    const a = cp[i], b = cp[i+1];
    return [a[0]+(b[0]-a[0])*u, a[1]+(b[1]-a[1])*u, a[2]+(b[2]-a[2])*u];
}

/* ---------- WebGL ---------- */
function compileShader(gl, src, type) {
    const s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
        console.error('Shader compile error:', gl.getShaderInfoLog(s));
        gl.deleteShader(s);
        return null;
    }
    return s;
}

function createProgram(gl, vs, fs) {
    const p = gl.createProgram();
    gl.attachShader(p, vs);
    gl.attachShader(p, fs);
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
        console.error('Link error:', gl.getProgramInfoLog(p));
        return null;
    }
    return p;
}

async function initWebGL(fsSrc) {
    const canvas = document.getElementById('canvasGL');
    gl = canvas.getContext('webgl2', {preserveDrawingBuffer:false});
    if (!gl) { alert('WebGL2 not supported'); return; }

    const vsSrc = document.getElementById('vs').text.trim();
    const vs = compileShader(gl, vsSrc, gl.VERTEX_SHADER);
    const fs = compileShader(gl, fsSrc.trim(), gl.FRAGMENT_SHADER);
    program = createProgram(gl, vs, fs);

    const buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1,-1,1,-1,-1,1,1,1]), gl.STATIC_DRAW);

    const posLoc = gl.getAttribLocation(program, 'a_position');
    gl.enableVertexAttribArray(posLoc);
    gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);

    locs = {
        resolution: gl.getUniformLocation(program, 'u_resolution'),
        frag1_pos:  gl.getUniformLocation(program, 'u_frag1_pos'),
        frag1_reqh: gl.getUniformLocation(program, 'u_frag1_reqh'),
        frag2_pos:  gl.getUniformLocation(program, 'u_frag2_pos'),
        frag2_reqh: gl.getUniformLocation(program, 'u_frag2_reqh'),
        n1: gl.getUniformLocation(program, 'u_n1'),
        n2: gl.getUniformLocation(program, 'u_n2'),
        pivot1: gl.getUniformLocation(program, 'u_pivot1'),
        angleRange: gl.getUniformLocation(program, 'u_angleRange'),
        distRange: gl.getUniformLocation(program, 'u_distRange'),
        mode: gl.getUniformLocation(program, 'u_mode'),
        cartRange: gl.getUniformLocation(program, 'u_cartRange'),
        mask: gl.getUniformLocation(program, 'u_mask'),
        b: gl.getUniformLocation(program, 'u_b'),
        colormap: gl.getUniformLocation(program, 'u_colormap'),
        energyMin: gl.getUniformLocation(program, 'u_energyMin'),
        energyMax: gl.getUniformLocation(program, 'u_energyMax')
    };

    buildColormaps();
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    updateColormap();

    resizeCanvas();
    window.addEventListener('resize', resizeCanvas);
}

let currentColormap = 'seismic';
let colormapTex;

function updateColormap() {
    if (!colormapTex) {
        colormapTex = gl.createTexture();
        gl.bindTexture(gl.TEXTURE_2D, colormapTex);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    }
    const data = COLORMAPS[currentColormap];
    gl.bindTexture(gl.TEXTURE_2D, colormapTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGB8, 256, 1, 0, gl.RGB, gl.UNSIGNED_BYTE, data);
}

function resizeCanvas() {
    const canvas = document.getElementById('canvasGL');
    const rect = canvas.parentElement.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvas.width = rect.width * dpr;
    canvas.height = rect.height * dpr;
    gl.viewport(0, 0, canvas.width, canvas.height);

    const overlay = document.getElementById('overlay');
    overlay.width = rect.width * dpr;
    overlay.height = rect.height * dpr;

    const lp = document.getElementById('lineProfile');
    if (lp) {
        lp.width = rect.width * dpr;
        lp.height = 150 * dpr;
    }

    render();
}

/* ---------- Tick helpers ---------- */
function niceTicks(min, max, targetCount) {
    const span = max - min;
    if (span <= 0) return [];
    const rough = span / targetCount;
    const pow10 = Math.pow(10, Math.floor(Math.log10(rough)));
    const mults = [1, 2, 5, 10];
    let step = pow10 * 10;
    for (const m of mults) {
        if (m * pow10 >= rough * 0.8) { step = m * pow10; break; }
    }
    const t0 = Math.ceil(min / step) * step;
    const ticks = [];
    for (let t = t0; t <= max + 1e-9; t += step) ticks.push(t);
    return ticks;
}

function fmt(v) {
    if (Math.abs(v) < 1e-6) return '0';
    if (Math.abs(v) >= 100) return Math.round(v).toString();
    if (Math.abs(v) >= 10) return v.toFixed(1);
    return v.toFixed(2);
}

/* ---------- 1D Line Profile ---------- */
function drawLineProfile(profile, emin, emax) {
    const canvas = document.getElementById('lineProfile');
    if (!canvas || !profile || profile.length < 2) return;
    const ctx = canvas.getContext('2d');
    const w = canvas.width;
    const h = canvas.height;
    ctx.clearRect(0, 0, w, h);

    const padL = 50, padR = 15, padT = 15, padB = 28;
    const graphW = w - padL - padR;
    const graphH = h - padT - padB;

    const dMin = profile[0].y;
    const dMax = profile[profile.length - 1].y;
    const eRange = emax - emin;
    if (eRange < 1e-6) return;

    const xForD = (d) => padL + (d - dMin) / (dMax - dMin) * graphW;
    const yForE = (e) => padT + graphH - (e - emin) / eRange * graphH;

    // axes
    ctx.strokeStyle = '#666';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padL, padT);
    ctx.lineTo(padL, h - padB);
    ctx.lineTo(w - padR, h - padB);
    ctx.stroke();

    // labels
    ctx.fillStyle = '#aaa';
    ctx.font = '10px Arial';
    ctx.textAlign = 'center';
    ctx.fillText('Distance (Å)', padL + graphW / 2, h - 4);
    ctx.save();
    ctx.translate(10, padT + graphH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText('E (kcal/mol)', 0, 0);
    ctx.restore();

    // x-ticks (distance)
    const xTicks = niceTicks(dMin, dMax, 5);
    ctx.fillStyle = '#aaa';
    ctx.font = '9px Arial';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (const t of xTicks) {
        const x = xForD(t);
        ctx.strokeStyle = '#666';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(x, h - padB);
        ctx.lineTo(x, h - padB + 4);
        ctx.stroke();
        ctx.fillText(fmt(t), x, h - padB + 6);
    }

    // y-ticks (energy)
    const yTicks = niceTicks(emin, emax, 4);
    ctx.fillStyle = '#aaa';
    ctx.font = '9px Arial';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const t of yTicks) {
        const y = yForE(t);
        ctx.strokeStyle = '#666';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(padL, y);
        ctx.lineTo(padL - 4, y);
        ctx.stroke();
        ctx.fillText(fmt(t), padL - 6, y);
    }

    // plot line
    ctx.beginPath();
    ctx.moveTo(xForD(profile[0].y), yForE(profile[0].E));
    for (let i = 1; i < profile.length; i++) {
        ctx.lineTo(xForD(profile[i].y), yForE(profile[i].E));
    }
    ctx.strokeStyle = '#4CAF50';
    ctx.lineWidth = 2;
    ctx.stroke();

    // zero-energy dashed line
    ctx.strokeStyle = '#666';
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    const y0 = yForE(0);
    ctx.beginPath();
    ctx.moveTo(padL, y0);
    ctx.lineTo(w - padR, y0);
    ctx.stroke();
    ctx.setLineDash([]);
}

/* ---------- Overlay (debug drawing + axes) ---------- */
function drawOverlay() {
    const overlay = document.getElementById('overlay');
    const ctx = overlay.getContext('2d');
    const w = overlay.width;
    const h = overlay.height;
    ctx.clearRect(0, 0, w, h);

    if (!current.frag1) return;
    const mode = document.querySelector('input[name="mode"]:checked').value;

    ctx.save();
    ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

    const displayW = w / window.devicePixelRatio;
    const displayH = h / window.devicePixelRatio;
    const cartRange = parseFloat(document.getElementById('cartRange').value);
    const scale = Math.min(displayW, displayH) / (2 * cartRange);
    const cx = displayW / 2;
    const cy = displayH / 2;

    const colors = {O:'#ff4444', H:'#ffffff', C:'#aaaaaa', N:'#4444ff', F:'#44ff44', Lp:'#ff44ff'};
    const radii  = {O:0.15, H:0.08, C:0.15, N:0.15, F:0.15, Lp:0.08};

    const proj = document.getElementById('projSelect').value;
    const permute = (x, y, z) => {
        if (proj === 'xy') return [x, y, z];
        if (proj === 'xz') return [x, z, -y];
        if (proj === 'yz') return [y, z, -x];
        if (proj === 'zy') return [z, y, -x];
        return [x, y, z];
    };

    if (mode === 'cartesian') {
        // Cartesian axes through center with world-coordinate ticks
        ctx.strokeStyle = '#666';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(0, cy); ctx.lineTo(displayW, cy);
        ctx.moveTo(cx, 0); ctx.lineTo(cx, displayH);
        ctx.stroke();

        const xTicks = niceTicks(-cartRange, cartRange, 6);
        ctx.fillStyle = '#aaa';
        ctx.font = '9px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        for (const t of xTicks) {
            const px = cx + t * scale;
            ctx.beginPath();
            ctx.moveTo(px, cy - 3);
            ctx.lineTo(px, cy + 3);
            ctx.strokeStyle = '#666';
            ctx.stroke();
            ctx.fillText(fmt(t), px, cy + 5);
        }
        const yTicks = niceTicks(-cartRange, cartRange, 6);
        ctx.fillStyle = '#aaa';
        ctx.font = '9px Arial';
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        for (const t of yTicks) {
            const py = cy - t * scale;
            ctx.beginPath();
            ctx.moveTo(cx - 3, py);
            ctx.lineTo(cx + 3, py);
            ctx.strokeStyle = '#666';
            ctx.stroke();
            ctx.fillText(fmt(t), cx - 5, py);
        }
    } else {
        // Polar mode: edge ticks for angle (bottom) and distance (left)
        const angleMinDeg = parseFloat(document.getElementById('angleMin').value);
        const angleMaxDeg = parseFloat(document.getElementById('angleMax').value);
        const dMin = parseFloat(document.getElementById('distMin').value);
        const dMax = parseFloat(document.getElementById('distMax').value);

        // Bottom edge: angle ticks
        const aTicksDeg = niceTicks(angleMinDeg, angleMaxDeg, 6);
        ctx.fillStyle = '#aaa';
        ctx.font = '9px Arial';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.strokeStyle = '#666';
        ctx.lineWidth = 1;
        for (const t of aTicksDeg) {
            const frac = (t - angleMinDeg) / (angleMaxDeg - angleMinDeg);
            const px = frac * displayW;
            ctx.beginPath();
            ctx.moveTo(px, displayH - 4);
            ctx.lineTo(px, displayH);
            ctx.stroke();
            ctx.fillText(Math.round(t) + '\u00B0', px, displayH - 8);
        }

        // Left edge: distance ticks
        const dTicks = niceTicks(dMin, dMax, 5);
        ctx.fillStyle = '#aaa';
        ctx.font = '9px Arial';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'middle';
        ctx.strokeStyle = '#666';
        for (const t of dTicks) {
            const frac = (t - dMin) / (dMax - dMin);
            const py = displayH - frac * displayH;
            ctx.beginPath();
            ctx.moveTo(0, py);
            ctx.lineTo(4, py);
            ctx.stroke();
            ctx.fillText(fmt(t), 6, py);
        }
    }

    const ppx = current.frag1.atoms[current.pivot1].x;
    const ppy = current.frag1.atoms[current.pivot1].y;
    const ppz = current.frag1.atoms[current.pivot1].z;

    if (mode === 'cartesian') {
        // Draw fragment 1 atoms in cartesian coordinates
        for (let i = 0; i < current.frag1.n; i++) {
            const a = current.frag1.atoms[i];
            const [pxRaw, pyRaw, pzRaw] = permute(a.x - ppx, a.y - ppy, a.z - ppz);
            const px = cx + pxRaw * scale;
            const py = cy - pyRaw * scale;
            const r = (radii[a.el] || 0.1) * scale;

            ctx.beginPath();
            ctx.arc(px, py, r, 0, 2*Math.PI);
            ctx.fillStyle = colors[a.el] || '#ccc';
            ctx.fill();
            ctx.strokeStyle = '#000';
            ctx.lineWidth = 1;
            ctx.stroke();

            ctx.fillStyle = '#fff';
            ctx.font = 'bold 11px Arial';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText(a.el + (i+1), px, py);
        }
    } else {
        // Draw fragment 1 atoms in polar coordinates (R, theta)
        const angleMinDeg = parseFloat(document.getElementById('angleMin').value);
        const angleMaxDeg = parseFloat(document.getElementById('angleMax').value);
        const aMin = angleMinDeg * Math.PI / 180;
        const aMax = angleMaxDeg * Math.PI / 180;
        const dMin = parseFloat(document.getElementById('distMin').value);
        const dMax = parseFloat(document.getElementById('distMax').value);

        for (let i = 0; i < current.frag1.n; i++) {
            const a = current.frag1.atoms[i];
            const [x, y, z] = permute(a.x - ppx, a.y - ppy, a.z - ppz);
            const R = Math.sqrt(x*x + y*y);
            const theta = Math.atan2(x, y); // angle from +Y toward +X, matching shader
            let px = displayW * (theta - aMin) / (aMax - aMin);
            let py = displayH * (dMax - R) / (dMax - dMin);

            // clip to plot area so every atom is visible
            const onScreen = (px >= -5 && px <= displayW + 5 && py >= -5 && py <= displayH + 5);
            px = Math.max(2, Math.min(displayW - 2, px));
            py = Math.max(2, Math.min(displayH - 2, py));

            // small cross marker
            const r = 3;
            ctx.strokeStyle = colors[a.el] || '#ccc';
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.moveTo(px - r, py - r); ctx.lineTo(px + r, py + r);
            ctx.moveTo(px - r, py + r); ctx.lineTo(px + r, py - r);
            ctx.stroke();

            // black label with tiny white halo for readability
            ctx.font = 'bold 10px Arial';
            ctx.textAlign = 'left';
            ctx.textBaseline = 'bottom';
            const label = a.el + (i+1);
            ctx.fillStyle = '#fff';
            ctx.fillText(label, px + 2, py - 1);
            ctx.fillStyle = '#000';
            ctx.fillText(label, px + 1, py - 2);
        }
    }

    if (mode === 'cartesian' && current.debugClick && current.frag2) {
        const dc = current.debugClick;
        const clickX = dc.x / window.devicePixelRatio;
        const clickY = dc.y / window.devicePixelRatio;
        const worldX = -cartRange + (2 * cartRange) * (clickX / displayW);
        const worldY = cartRange - (2 * cartRange) * (clickY / displayH);
        const toPixel = [worldX, worldY, 0];
        const tpLen = Math.sqrt(toPixel[0]*toPixel[0] + toPixel[1]*toPixel[1]);
        let dir, perp;
        if (tpLen < 0.001) {
            dir = [1, 0]; perp = [0, 1];
        } else {
            dir = [toPixel[0]/tpLen, toPixel[1]/tpLen];
            perp = [-dir[1], dir[0]];
        }
        for (let j = 0; j < current.frag2.n; j++) {
            const a = current.frag2.atoms[j];
            const px = cx + (worldX + dir[0]*a.x + perp[0]*a.y) * scale;
            const py = cy - (worldY + dir[1]*a.x + perp[1]*a.y) * scale;
            const r = (radii[a.el] || 0.1) * scale;
            ctx.beginPath();
            ctx.arc(px, py, r, 0, 2*Math.PI);
            ctx.fillStyle = colors[a.el] || '#ccc';
            ctx.globalAlpha = 0.7;
            ctx.fill();
            ctx.globalAlpha = 1.0;
            ctx.strokeStyle = '#ff0';
            ctx.lineWidth = 1;
            ctx.stroke();
            ctx.fillStyle = '#ff0';
            ctx.font = 'bold 10px Arial';
            ctx.fillText(a.el + (j+1), px, py);
        }
    }

    ctx.restore();
}

/* ---------- GUI ---------- */
function populateSelect(id, names) {
    const sel = document.getElementById(id);
    sel.innerHTML = '';
    names.forEach(n => {
        const o = document.createElement('option');
        o.value = n; o.textContent = n;
        sel.appendChild(o);
    });
}

function makeTable(containerId, fragObj, prefix) {
    const div = document.getElementById(containerId);
    div.innerHTML = '';
    const tbl = document.createElement('table');
    tbl.className = 'req-table';
    const hdr = document.createElement('tr');
    ['Atom','R','E','Q','H'].forEach(h => {
        const th = document.createElement('th');
        th.textContent = h;
        hdr.appendChild(th);
    });
    tbl.appendChild(hdr);

    fragObj.atoms.forEach((a, i) => {
        const tr = document.createElement('tr');
        const tdLabel = document.createElement('td');
        tdLabel.className = 'atom-label';
        tdLabel.textContent = a.el + (i+1);
        tr.appendChild(tdLabel);

        ['R','E','Q','H'].forEach(key => {
            const td = document.createElement('td');
            const inp = document.createElement('input');
            inp.type = 'number';
            inp.value = a[key].toFixed(key==='E'||key==='Q'?3:2);
            inp.step = key==='E'||key==='Q'?0.001:0.01;
            inp.dataset.frag = prefix;
            inp.dataset.idx = i;
            inp.dataset.key = key;
            inp.addEventListener('change', () => {
                a[key] = parseFloat(inp.value) || 0;
                render();
            });
            inp.addEventListener('wheel', (e) => {
                e.preventDefault();
                const step = parseFloat(inp.step) || 0.01;
                a[key] += (e.deltaY < 0 ? 1 : -1) * step;
                inp.value = a[key].toFixed(key==='E'||key==='Q'?3:2);
                render();
            });
            td.appendChild(inp);
            tr.appendChild(td);
        });
        tbl.appendChild(tr);
    });
    div.appendChild(tbl);
}

function updatePivotSelect() {
    const sel = document.getElementById('pivotSelect');
    sel.innerHTML = '';
    if (!current.frag1) return;
    current.frag1.atoms.forEach((a, i) => {
        const o = document.createElement('option');
        o.value = i;
        o.textContent = a.el + (i+1);
        if (i === current.pivot1) o.selected = true;
        sel.appendChild(o);
    });
}

async function switchFrag1(name) {
    current.frag1 = await loadMolecule(name);
    current.pivot1 = 0;
    updatePivotSelect();
    makeTable('frag1Table', current.frag1, 'f1');
    render();
}

async function switchFrag2(name) {
    current.frag2 = await loadMolecule(name);
    makeTable('frag2Table', current.frag2, 'f2');
    render();
}

/* ---------- Render ---------- */
function render() {
    if (!gl || !program || !current.frag1 || !current.frag2) return;

    const canvas = document.getElementById('canvasGL');
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.useProgram(program);

    const pivotSel = document.getElementById('pivotSelect');
    current.pivot1 = parseInt(pivotSel.value) || 0;
    const pivotAtom = current.frag1.atoms[current.pivot1];
    const ppx = pivotAtom.x, ppy = pivotAtom.y, ppz = pivotAtom.z;

    const proj = document.getElementById('projSelect').value;
    const permute = (x, y, z) => {
        if (proj === 'xy') return [x, y, z];
        if (proj === 'xz') return [x, z, -y];
        if (proj === 'yz') return [y, z, -x];
        if (proj === 'zy') return [z, y, -x];
        return [x, y, z];
    };

    const f1p = new Float32Array(10*4);
    const f1r = new Float32Array(10*4);
    current.frag1.atoms.forEach((a, i) => {
        const [px, py, pz] = permute(a.x - ppx, a.y - ppy, a.z - ppz);
        f1p[i*4] = px; f1p[i*4+1] = py; f1p[i*4+2] = pz; f1p[i*4+3] = 0;
        f1r[i*4] = a.R; f1r[i*4+1] = a.E; f1r[i*4+2] = a.Q; f1r[i*4+3] = a.H;
    });
    const f2p = new Float32Array(10*4);
    const f2r = new Float32Array(10*4);
    current.frag2.atoms.forEach((a, i) => {
        f2p[i*4] = a.x; f2p[i*4+1] = a.y; f2p[i*4+2] = a.z; f2p[i*4+3] = 0;
        f2r[i*4] = a.R; f2r[i*4+1] = a.E; f2r[i*4+2] = a.Q; f2r[i*4+3] = a.H;
    });

    gl.uniform4fv(locs.frag1_pos, f1p);
    gl.uniform4fv(locs.frag1_reqh, f1r);
    gl.uniform4fv(locs.frag2_pos, f2p);
    gl.uniform4fv(locs.frag2_reqh, f2r);
    gl.uniform1i(locs.n1, current.frag1.n);
    gl.uniform1i(locs.n2, current.frag2.n);
    gl.uniform1i(locs.pivot1, current.pivot1);

    const mode = document.querySelector('input[name="mode"]:checked').value;
    gl.uniform1i(locs.mode, mode === 'cartesian' ? 1 : 0);
    const cartRange = parseFloat(document.getElementById('cartRange').value);
    gl.uniform1f(locs.cartRange, cartRange);

    const aMin = parseFloat(document.getElementById('angleMin').value) * Math.PI / 180;
    const aMax = parseFloat(document.getElementById('angleMax').value) * Math.PI / 180;
    const dMin = parseFloat(document.getElementById('distMin').value);
    const dMax = parseFloat(document.getElementById('distMax').value);
    gl.uniform2f(locs.angleRange, aMin, aMax);
    gl.uniform2f(locs.distRange, dMin, dMax);

    const mask = [
        document.getElementById('maskPauli').checked ? 1 : 0,
        document.getElementById('maskLondon').checked ? 1 : 0,
        document.getElementById('maskCoulomb').checked ? 1 : 0,
        document.getElementById('maskHB').checked ? 1 : 0
    ];
    gl.uniform4f(locs.mask, mask[0], mask[1], mask[2], mask[3]);

    const b = parseFloat(document.getElementById('bParam').value);
    gl.uniform1f(locs.b, b);

    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, colormapTex);
    gl.uniform1i(locs.colormap, 0);

    const auto = document.getElementById('autoRange').checked;
    let vmax = parseFloat(document.getElementById('vmax').value);
    const profile = Molecule.lineProfile(current.frag1, current.frag2, current.pivot1, mode, aMin, aMax, dMin, dMax, cartRange, b, mask, 64);
    if (auto) {
        let vmin = 1e9;
        for (const p of profile) { if (p.E < vmin) vmin = p.E; }
        vmax = Math.round(Math.abs(vmin));
        if (vmax < 0.1) vmax = 0.1;
        document.getElementById('vmax').value = vmax;
    }
    if (vmax < 0.1) vmax = 0.1;
    document.getElementById('vmaxVal').textContent = vmax;
    gl.uniform1f(locs.energyMin, -vmax);
    gl.uniform1f(locs.energyMax, vmax);

    gl.uniform2f(locs.resolution, canvas.width, canvas.height);

    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    drawOverlay();
    drawLineProfile(profile, -vmax, vmax);
}

/* ---------- Init ---------- */
async function init() {
    let fsSrc;
    try {
        const r = await fetch('shader.frag');
        fsSrc = await r.text();
    } catch(e) {
        alert('Cannot load shader.frag. Please serve this folder via HTTP (e.g. python -m http.server).');
        return;
    }

    await initWebGL(fsSrc);

    populateSelect('frag1Select', Object.keys(MOLECULES).filter(k => MOLECULES[k].type==='frag1'));
    populateSelect('frag2Select', Object.keys(MOLECULES).filter(k => MOLECULES[k].type==='frag2'));

    const cmSel = document.getElementById('colormapSelect');
    Object.keys(COLORMAPS).forEach(name => {
        const o = document.createElement('option');
        o.value = name; o.textContent = name;
        if (name === 'seismic') o.selected = true;
        cmSel.appendChild(o);
    });

    document.getElementById('frag1Select').addEventListener('change', e => switchFrag1(e.target.value));
    document.getElementById('frag2Select').addEventListener('change', e => switchFrag2(e.target.value));
    document.getElementById('pivotSelect').addEventListener('change', () => render());
    document.getElementById('projSelect').addEventListener('change', () => render());
    ['angleMin','angleMax','distMin','distMax','bParam','vmax','cartRange'].forEach(id => {
        document.getElementById(id).addEventListener('input', () => {
            if (id==='bParam') document.getElementById('bVal').textContent = document.getElementById(id).value;
            if (id==='vmax') document.getElementById('vmaxVal').textContent = document.getElementById(id).value;
            if (id==='cartRange') document.getElementById('cartVal').textContent = document.getElementById(id).value;
            render();
        });
    });
    document.querySelectorAll('input[name="mode"]').forEach(r => {
        r.addEventListener('change', () => render());
    });
    ['maskPauli','maskLondon','maskCoulomb','maskHB','autoRange','addEpair'].forEach(id => {
        document.getElementById(id).addEventListener('change', () => {
            if (id === 'addEpair') {
                switchFrag1(document.getElementById('frag1Select').value);
            } else {
                render();
            }
        });
    });
    document.getElementById('lpDist').addEventListener('input', () => {
        document.getElementById('lpDistVal').textContent = document.getElementById('lpDist').value;
        if (document.getElementById('addEpair').checked) {
            switchFrag1(document.getElementById('frag1Select').value);
        }
    });
    cmSel.addEventListener('change', e => {
        currentColormap = e.target.value;
        updateColormap();
        render();
    });

    const overlay = document.getElementById('overlay');
    overlay.addEventListener('click', (e) => {
        const rect = overlay.getBoundingClientRect();
        const x = (e.clientX - rect.left) * window.devicePixelRatio;
        const y = (e.clientY - rect.top) * window.devicePixelRatio;
        current.debugClick = {x, y, rectW: rect.width * window.devicePixelRatio, rectH: rect.height * window.devicePixelRatio};
        render();
    });

    document.getElementById('btnQEq1').addEventListener('click', () => {
        if (!current.frag1) return;
        runQEqOnFragment(current.frag1.atoms, 0.0);
        makeTable('frag1Table', current.frag1, 'f1');
        render();
    });
    document.getElementById('btnQEq2').addEventListener('click', () => {
        if (!current.frag2) return;
        runQEqOnFragment(current.frag2.atoms, 0.0);
        makeTable('frag2Table', current.frag2, 'f2');
        render();
    });

    // Default to CH2NH with pivot on N (atom 1)
    document.getElementById('frag1Select').value = 'CH2NH';
    await switchFrag1('CH2NH');
    current.pivot1 = 1;
    updatePivotSelect();

    await switchFrag2(document.getElementById('frag2Select').value);
    render();
}

init();
