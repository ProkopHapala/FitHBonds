#version 300 es
precision highp float;

in vec2 v_uv;
out vec4 fragColor;

uniform vec2  u_resolution;
uniform vec4  u_frag1_pos[10];
uniform vec4  u_frag1_reqh[10];
uniform vec4  u_frag2_pos[10];
uniform vec4  u_frag2_reqh[10];
uniform int   u_n1;
uniform int   u_n2;
uniform int   u_pivot1;
uniform vec2  u_angleRange;
uniform vec2  u_distRange;
uniform int   u_mode;       // 0=polar, 1=cartesian
uniform float u_cartRange;
uniform vec4  u_mask;
uniform float u_b;
uniform sampler2D u_colormap;
uniform float u_energyMin;
uniform float u_energyMax;

float hbCorrection(int i, int j, float r, vec4 ri, vec4 rj) {
    return 0.0;
}

void main() {
    vec3 center;
    vec3 dir;
    vec3 perp;
    vec3 normal = vec3(0.0, 0.0, 1.0);
    vec3 pivotPos = u_frag1_pos[u_pivot1].xyz;

    if (u_mode == 0) {
        // Polar mode: x=angle, y=distance
        float angle = mix(u_angleRange.x, u_angleRange.y, v_uv.x);
        float dist  = mix(u_distRange.x,  u_distRange.y,  v_uv.y);
        float ca = cos(angle);
        float sa = sin(angle);
        dir   = vec3(sa, ca, 0.0);
        perp  = vec3(-ca, sa, 0.0);
        center = pivotPos + dist * dir;
    } else {
        // Cartesian mode: pixel position is world position
        float worldX = mix(-u_cartRange, u_cartRange, v_uv.x);
        float worldY = mix(-u_cartRange, u_cartRange, v_uv.y);
        center = vec3(worldX, worldY, 0.0);
        vec3 toPixel = center - pivotPos;
        float tpLen = length(toPixel);
        if (tpLen < 0.001) {
            dir = vec3(1.0, 0.0, 0.0);
        } else {
            dir = toPixel / tpLen;
        }
        perp = vec3(-dir.y, dir.x, 0.0);
    }

    float totalEnergy = 0.0;

    for (int i = 0; i < 10; i++) {
        if (i >= u_n1) break;
        vec4 pi = u_frag1_pos[i];
        vec4 ri = u_frag1_reqh[i];
        float Ri = ri.x;
        float Ei = ri.y;
        float Qi = ri.z;

        for (int j = 0; j < 10; j++) {
            if (j >= u_n2) break;
            vec4 pj_local = u_frag2_pos[j];
            vec4 rj = u_frag2_reqh[j];
            float Rj = rj.x;
            float Ej = rj.y;
            float Qj = rj.z;

            vec3 pj_world = center
                + dir   * pj_local.x
                + perp  * pj_local.y
                + normal * pj_local.z;

            vec3 d = pi.xyz - pj_world;
            float r = length(d);
            if (r < 0.05) r = 0.05;

            float Rij = Ri + Rj;
            float Eij = sqrt(max(0.0, Ei * Ej));
            float dr  = r - Rij;
            float eb  = exp(-u_b * dr);
            float eb2 = eb * eb;

            if (u_mask.x > 0.5) {
                totalEnergy += Eij * eb2;
            }
            if (u_mask.y > 0.5) {
                totalEnergy -= 2.0 * Eij * eb;
            }
            if (u_mask.z > 0.5) {
                totalEnergy += 332.0637 * Qi * Qj / r;
            }
            if (u_mask.w > 0.5) {
                totalEnergy += hbCorrection(i, j, r, ri, rj);
            }
        }
    }

    float t = (totalEnergy - u_energyMin) / (u_energyMax - u_energyMin);
    t = clamp(t, 0.0, 1.0);

    vec3 color = texture(u_colormap, vec2(t, 0.5)).rgb;
    fragColor = vec4(color, 1.0);
}
