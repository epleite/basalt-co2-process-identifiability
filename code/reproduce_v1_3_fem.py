from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'frozen'
OUT = ROOT / 'reproduced' / 'v1_3'
OUT.mkdir(parents=True, exist_ok=True)

MAT = pd.read_csv(DATA / 'basalt_co2_v1_1_nominal_SIP_parameters.csv')
MAT = MAT[MAT.model_form == 'DEM-tangent'].copy()

FREQS = np.array([0.10, 0.30, 1.0, 3.0])
TIMES = [365, 1825]
FLOW = {'upper': 9.507675310035554, 'core': 7.573390737497025, 'fracture': 1.7338595036951663}
H = sum(FLOW.values())
WELL_X = [0.0, 25.0]
Z_E = np.unique(np.r_[np.arange(-40, -20, 5.0), np.arange(-20, 20.1, 2.5), np.arange(25, 40.1, 5.0)])
BG_RHO, BG_MN, BG_TAU, BG_C = 150.0, 5e-4, .8, .4


def sigma_from_row(row, f):
    sig0 = 1.0 / row.rhoDC_ohm_m
    Mn, tau, c = row.Mn_S_per_m, row.tau_s, row.cole_c
    siginf = sig0 + Mn
    M = Mn / siginf
    return siginf * (1 - M / (1 + (1j * 2 * np.pi * f * tau) ** c))


def sigma_mat(t, proc, fac, f):
    r = MAT[(MAT.time_days == t) & (MAT.process == proc) & (MAT.facies == fac)].iloc[0]
    return sigma_from_row(r, f)


def facies_at(z):
    top = -H / 2
    if top <= z < top + FLOW['upper']:
        return 'upper'
    if top + FLOW['upper'] <= z < top + FLOW['upper'] + FLOW['core']:
        return 'core'
    if top + FLOW['upper'] + FLOW['core'] <= z <= H / 2:
        return 'fracture'
    return None


def sigma_bg(f):
    s0 = 1 / BG_RHO
    si = s0 + BG_MN
    M = BG_MN / si
    return si * (1 - M / (1 + (1j * 2 * np.pi * f * BG_TAU) ** BG_C))


# Mesh
DX = DZ = 1.25
XS = np.arange(-15.0, 40.0 + 1e-9, DX)
ZS = np.arange(-50.0, 50.0 + 1e-9, DZ)
NX, NZ = len(XS), len(ZS)
X, Z = np.meshgrid(XS, ZS)
NODES = np.c_[X.ravel(), Z.ravel()]

def nid(ix, iz): return iz * NX + ix
TRIS = []
for iz in range(NZ - 1):
    for ix in range(NX - 1):
        n00, n10, n01, n11 = nid(ix, iz), nid(ix + 1, iz), nid(ix, iz + 1), nid(ix + 1, iz + 1)
        if (ix + iz) % 2 == 0:
            TRIS += [(n00, n10, n11), (n00, n11, n01)]
        else:
            TRIS += [(n00, n10, n01), (n10, n11, n01)]
TRIS = np.asarray(TRIS, int)
NE = len(TRIS)
AREAS = np.empty(NE)
BGRAD = np.empty((NE, 2, 3))
CENT = np.empty((NE, 2))
for e, t in enumerate(TRIS):
    xy = NODES[t]
    x1, z1 = xy[0]; x2, z2 = xy[1]; x3, z3 = xy[2]
    A = .5 * ((x2 - x1) * (z3 - z1) - (x3 - x1) * (z2 - z1))
    if A < 0:
        t = t[[0, 2, 1]]; TRIS[e] = t; xy = NODES[t]
        x1, z1 = xy[0]; x2, z2 = xy[1]; x3, z3 = xy[2]
        A = .5 * ((x2 - x1) * (z3 - z1) - (x3 - x1) * (z2 - z1))
    AREAS[e] = A
    b = np.array([z2 - z3, z3 - z1, z1 - z2]) / (2 * A)
    c = np.array([x3 - x2, x1 - x3, x2 - x1]) / (2 * A)
    BGRAD[e] = np.vstack([b, c])
    CENT[e] = xy.mean(axis=0)
KE_GEOM = np.asarray([AREAS[e] * (BGRAD[e].T @ BGRAD[e]) for e in range(NE)])
BOUNDARY = (np.isclose(NODES[:, 0], XS[0]) | np.isclose(NODES[:, 0], XS[-1]) |
            np.isclose(NODES[:, 1], ZS[0]) | np.isclose(NODES[:, 1], ZS[-1]))
FREE = np.where(~BOUNDARY)[0]
FREE_MAP = -np.ones(len(NODES), int); FREE_MAP[FREE] = np.arange(len(FREE))


def electrode_group(x, z0):
    ix = np.argmin(np.abs(XS - x)); iz0 = np.argmin(np.abs(ZS - z0))
    return np.asarray([nid(ix, iz0 + off) for off in [-1, 0, 1] if 0 <= iz0 + off < NZ], int)

ELECS = []
for wi, x in enumerate(WELL_X):
    for z0 in Z_E:
        ELECS.append((wi, float(z0), electrode_group(x, z0)))
NELEC = len(ELECS)


def electrode_rhs(eid):
    rhs = np.zeros(len(FREE), complex)
    grp = ELECS[eid][2]
    valid = [n for n in grp if FREE_MAP[n] >= 0]
    for n in valid:
        rhs[FREE_MAP[n]] += 1 / len(valid)
    return rhs
BSRC = np.column_stack([electrode_rhs(i) for i in range(NELEC)])


def assemble(sig):
    rows, cols, vals = [], [], []
    for e, t in enumerate(TRIS):
        Ke = sig[e] * KE_GEOM[e]
        for a in range(3):
            ia = FREE_MAP[t[a]]
            if ia < 0: continue
            for b in range(3):
                ib = FREE_MAP[t[b]]
                if ib < 0: continue
                rows.append(ia); cols.append(ib); vals.append(Ke[a, b])
    return csc_matrix((vals, (rows, cols)), shape=(len(FREE), len(FREE)))


def solve_all(sig):
    lu = splu(assemble(sig))
    Ufree = lu.solve(BSRC)
    U = np.zeros((len(NODES), NELEC), complex); U[FREE] = Ufree
    Ve = np.zeros((NELEC, NELEC), complex)
    for j in range(NELEC):
        for i in range(NELEC):
            Ve[i, j] = U[ELECS[i][2], j].mean()
    return Ve


def element_sigma(t, proc, f, completion='resistive', plume_width=25.0):
    sig = np.full(NE, sigma_bg(f), dtype=complex)
    for e, (x, z) in enumerate(CENT):
        fac = facies_at(z)
        if fac is not None:
            sig[e] = sigma_mat(t, 'baseline', fac, f)
            if 0 <= x <= plume_width and proc != 'baseline':
                sig[e] = sigma_mat(t, proc, fac, f)
    if completion != 'ideal':
        scomp = (1 / (500.0 if completion == 'resistive' else 10.0)) + 0j
        for wx in WELL_X:
            sig[np.abs(CENT[:, 0] - wx) <= .75] = scomp
    return sig


def candidate_quads():
    left = list(range(len(Z_E))); right = list(range(len(Z_E), 2 * len(Z_E)))
    q = []
    for i in range(len(Z_E)):
        for j in range(i + 1, len(Z_E)):
            if 2.5 <= abs(Z_E[j] - Z_E[i]) <= 30:
                q.append((left[i], right[i], left[j], right[j]))
    lookup = {round(float(v), 3): i for i, v in enumerate(Z_E)}
    for sep in [5., 10., 20.]:
        pairs = []
        for i, z0 in enumerate(Z_E):
            key = round(float(z0 + sep), 3)
            if key in lookup: pairs.append((i, lookup[key]))
        for i, j in pairs:
            for k, l in pairs:
                q.append((left[i], left[j], right[k], right[l]))
    return q


def quad_data(Ve, qlist):
    return np.asarray([(Ve[M, A] - Ve[M, B]) - (Ve[N, A] - Ve[N, B]) for A, B, M, N in qlist])


def signed_phase(z):
    s = np.where(z.real >= 0, 1.0, -1.0)
    return np.angle(z * s)


def select_quads(all_quads):
    # FEM-baseline 1% target perturbation; independent of process identity
    scores = []
    baseline_by_f = []
    sens_by_f = []
    for f in FREQS:
        sig0 = element_sigma(1825, 'baseline', f, 'resistive', 25.)
        Ve0 = solve_all(sig0)
        db = quad_data(Ve0, all_quads)
        baseline_by_f.append(db)
        sigp = sig0.copy()
        target = (CENT[:, 0] >= 0) & (CENT[:, 0] <= 25) & (np.abs(CENT[:, 1]) <= H / 2)
        sigp[target] *= 1.01
        dp = quad_data(solve_all(sigp), all_quads)
        ph0, php = signed_phase(db), signed_phase(dp)
        dph = (php - ph0 + np.pi) % (2 * np.pi) - np.pi
        dlnamp = np.log(np.abs(dp / db))
        stable = (np.sign(db.real) == np.sign(dp.real)) & (np.abs(ph0) < .35) & (np.abs(php) < .50)
        score = np.sqrt((dph / .001) ** 2 + (dlnamp / .01) ** 2)
        score[~stable] = 0
        sens_by_f.append(score)
    score = np.sqrt(np.mean(np.asarray(sens_by_f) ** 2, axis=0))
    amp_stable = np.ones(len(all_quads), bool)
    for db in baseline_by_f:
        amp_stable &= np.abs(db) >= np.quantile(np.abs(db), .25)
    score[~amp_stable] = 0
    return [all_quads[i] for i in np.argsort(score)[-300:]]


def process_phase(t, proc, qlist):
    arr, stable = [], np.ones(len(qlist), bool)
    for f in FREQS:
        db = quad_data(solve_all(element_sigma(t, 'baseline', f, 'resistive', 25.)), qlist)
        dm = quad_data(solve_all(element_sigma(t, proc, f, 'resistive', 25.)), qlist)
        st = (np.sign(db.real) == np.sign(dm.real)) & (np.abs(signed_phase(db)) < .35) & (np.abs(signed_phase(dm)) < .50)
        stable &= st
        dph = (signed_phase(dm) - signed_phase(db) + np.pi) % (2 * np.pi) - np.pi
        arr.append(dph)
    return np.asarray(arr), stable


def main():
    allq = candidate_quads()
    q = select_quads(allq)
    rows = []
    cache = {}
    for t in TIMES:
        pp, stable = {}, np.ones(len(q), bool)
        for proc in ['dissolution', 'cementation', 'full']:
            arr, st = process_phase(t, proc, q)
            pp[proc] = arr
            stable &= st
        d = pp['dissolution'][:, stable].ravel()
        c = pp['cementation'][:, stable].ravel()
        cs = np.dot(d, c) / (np.linalg.norm(d) * np.linalg.norm(c))
        rmsD = np.sqrt(np.mean(d*d)) * 1000
        rmsC = np.sqrt(np.mean(c*c)) * 1000
        rmsF = np.sqrt(np.mean(pp['full'][:, stable].ravel() ** 2)) * 1000
        for pn in [1., 3., 5., 10.]:
            sd, sc = rmsD/pn, rmsC/pn
            s = min(sd, sc); w = s / np.sqrt(1 + s*s)
            seis = {365: -.9672, 1825: -.9974}[t]
            dc = {365: -.8623, 1825: -.9528}[t]
            joint = (seis + dc + w*w*cs) / (2 + w*w)
            rows.append([t, pn, int(stable.sum()), cs, abs(cs), rmsD, rmsC, rmsF, sd, sc, w, joint, abs(joint)])
    out = pd.DataFrame(rows, columns=['time_days','phase_noise_mrad','n_stable_quads','SIP_signed_cos','SIP_abs_cos',
                                     'phase_RMS_D_mrad','phase_RMS_C_mrad','phase_RMS_full_mrad','SNR_D','SNR_C',
                                     'SIP_weight','joint_signed_cos','joint_abs_cos'])
    out.to_csv(OUT / 'FEM_target_optimized_reproduced.csv', index=False)
    print(out.to_string(index=False))
    frozen = pd.read_csv(DATA / 'basalt_co2_v1_3_FEM_target_optimized.csv')
    merged = out.merge(frozen, on=['time_days','phase_noise_mrad'], suffixes=('_reproduced','_frozen'))
    cols = ['SIP_signed_cos','phase_RMS_D_mrad','phase_RMS_C_mrad','phase_RMS_full_mrad','joint_abs_cos']
    diffs = []
    for c in cols:
        diffs.append(np.max(np.abs(merged[c+'_reproduced'] - merged[c+'_frozen'])))
    print('max difference across headline FEM metrics:', max(diffs))


if __name__ == '__main__':
    main()
