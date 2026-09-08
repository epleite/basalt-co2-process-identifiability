from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data' / 'frozen'
OUT = ROOT / 'reproduced' / 'v1_4'
OUT.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(1414)


def gram_vectors(norms, gram):
    gram = np.asarray(gram, float)
    w, V = np.linalg.eigh(gram)
    w = np.clip(w, 1e-10, None)
    G = V @ np.diag(w) @ V.T
    L = np.linalg.cholesky(G)
    return np.asarray(norms)[:, None] * L


def canonical_3(c12, c13, c23, norms):
    G = np.array([[1, c12, c13], [c12, 1, c23], [c13, c23, 1]], float)
    return gram_vectors(norms, G)


def canonical_2(c, n1, n2, dim=3):
    v1 = np.zeros(dim)
    v2 = np.zeros(dim)
    v1[0] = n1
    v2[0] = n2 * c
    v2[1] = n2 * np.sqrt(max(1 - c * c, 0))
    return v1, v2


def pad(v, n):
    z = np.zeros(n)
    z[:len(v)] = v
    return z


def block(BS, BD, BC, Fref, name, time, modality):
    n = max(len(BS), len(BD), len(BC), len(Fref))
    BS, BD, BC, Fref = [pad(np.asarray(x, float), n) for x in [BS, BD, BC, Fref]]
    lin = BS + BD + BC
    return dict(B=np.column_stack([BS, BD, BC]), R=Fref - lin,
                name=name, time=time, modality=modality, full=Fref, lin=lin)


def weighted_quantile(x, w, qs):
    ix = np.argsort(x)
    xs = x[ix]
    ws = w[ix]
    c = np.cumsum(ws)
    c /= c[-1]
    return np.interp(qs, c, xs)


def weighted_stats(samples, w):
    w = w / w.sum()
    mean = np.sum(samples * w[:, None], axis=0)
    d = samples - mean
    cov = (d * w[:, None]).T @ d
    std = np.sqrt(np.diag(cov))
    corr = cov / np.outer(std, std)
    qs = np.array([weighted_quantile(samples[:, j], w, [.05, .5, .95])
                   for j in range(samples.shape[1])])
    return mean, std, cov, corr, qs[:, 0], qs[:, 1], qs[:, 2]


def normalize_logw(lw):
    m = np.max(lw)
    w = np.exp(lw - m)
    return w / w.sum()


def log_weights(pred, obs):
    r = pred - obs[None, :]
    return -0.5 * np.sum(r * r, axis=1)


def sigma_from_row(row, f):
    sig0 = 1.0 / row.rhoDC_ohm_m
    Mn = row.Mn_S_per_m
    tau = row.tau_s
    c = row.cole_c
    siginf = sig0 + Mn
    M = Mn / siginf
    return siginf * (1 - M / (1 + (1j * 2 * np.pi * f * tau) ** c))


def phase_vec_material(mat, t, proc, freqs):
    out = []
    for fac in ['upper', 'core', 'fracture']:
        rb = mat[(mat.time_days == t) & (mat.process == 'baseline') & (mat.facies == fac)].iloc[0]
        rp = mat[(mat.time_days == t) & (mat.process == proc) & (mat.facies == fac)].iloc[0]
        for f in freqs:
            out.append(np.angle(1 / sigma_from_row(rp, f)) - np.angle(1 / sigma_from_row(rb, f)))
    return np.asarray(out)


def cosine(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def build_blocks(phase_noise=5.0):
    v7 = pd.read_csv(DATA / 'basalt_co2_v0_7_basis_similarity.csv')
    e6 = pd.read_csv(DATA / 'basalt_co2_v0_6_reflectivity_energy.csv')
    dc_nom = pd.read_csv(DATA / 'basalt_co2_v1_0_nominal_resistivity.csv')
    fem = pd.read_csv(DATA / 'basalt_co2_v1_3_FEM_target_optimized.csv')
    procaware = pd.read_csv(DATA / 'basalt_co2_v1_3_process_aware_summary.csv')
    mat = pd.read_csv(DATA / 'basalt_co2_v1_1_nominal_SIP_parameters.csv')
    mat = mat[mat.model_form == 'DEM-tangent']

    blocks = {}
    # Seismic
    for t in [365, 1825]:
        z = v7[(v7.model_form == 'DEM-tangent') & (v7.time_days == t) &
               (v7.channel == 'PP+PS normalized')]
        cSD = float(z[(z.process_A == 'fluid') & (z.process_B == 'dissolution')].signed_cosine_similarity.iloc[0])
        cSC = float(z[(z.process_A == 'fluid') & (z.process_B == 'cementation')].signed_cosine_similarity.iloc[0])
        cDC = float(z[(z.process_A == 'dissolution') & (z.process_B == 'cementation')].signed_cosine_similarity.iloc[0])
        ee = e6[(e6.model_form == 'DEM-tangent') & (e6.wave_solver == 'exact-multiples') &
                (e6.acquisition == 'Geophone VSP') & (e6.time_days == t)].iloc[0]
        raw = np.sqrt([ee.E_fluid, ee.E_dissolution, ee.E_cementation])
        sigma = raw[2] / 4.0
        norms = raw / sigma
        V = canonical_3(cSD, cSC, cDC, norms)
        BS, BD, BC = V[0], V[1], V[2]
        lin = BS + BD + BC
        fullnorm = np.sqrt(ee.E_full) / sigma
        Fref = lin * (fullnorm / max(np.linalg.norm(lin), 1e-30))
        blocks[('seismic', t)] = block(BS, BD, BC, Fref, f'seis_{t}', t, 'seismic')

    # DC
    sigma_dc = .05
    for t in [365, 1825]:
        zz = dc_nom[(dc_nom.model_form == 'DEM-tangent') & (dc_nom.time_days == t)]
        def logvec(proc):
            r = zz[zz.process == proc].iloc[0]
            return np.log(np.array([r.rho_upper/r.rho0_upper,
                                    r.rho_core/r.rho0_core,
                                    r.rho_fracture/r.rho0_fracture])) / sigma_dc
        BS, BD, BC = [logvec(p) for p in ['fluid', 'dissolution', 'cementation']]
        blocks[('DC', t)] = block(BS, BD, BC, logvec('full'), f'dc_{t}', t, 'DC')

    # SIP
    freqs = np.array([.1, .3, 1., 3.])
    for t in [365, 1825]:
        fr = fem[(fem.time_days == t) & (fem.phase_noise_mrad == 5)].iloc[0]
        nD = fr.phase_RMS_D_mrad / phase_noise
        nC = fr.phase_RMS_C_mrad / phase_noise
        cDC = fr.SIP_signed_cos
        vm = {p: phase_vec_material(mat, t, p, freqs) for p in ['fluid', 'dissolution', 'cementation']}
        cSD = cosine(vm['fluid'], vm['dissolution'])
        cSC = cosine(vm['fluid'], vm['cementation'])
        if t == 365:
            scaleD = fr.phase_RMS_D_mrad / (np.sqrt(np.mean(vm['dissolution']**2)) * 1000)
            scaleC = fr.phase_RMS_C_mrad / (np.sqrt(np.mean(vm['cementation']**2)) * 1000)
            scale = .5 * (scaleD + scaleC)
            nS = (np.sqrt(np.mean(vm['fluid']**2)) * 1000 * scale) / phase_noise
            V = canonical_3(cSD, cSC, cDC, [nS, nD, nC])
            BS, BD, BC = [pad(v, 4) for v in V]
            lin = BS + BD + BC
            fullnorm = fr.phase_RMS_full_mrad / phase_noise
            Fref = lin.copy()
            Fref[-1] += np.sqrt(max(fullnorm**2 - np.linalg.norm(lin)**2, 0))
        else:
            BS = np.zeros(4)
            BD0, BC0 = canonical_2(cDC, nD, nC, dim=3)
            BD, BC = pad(BD0, 4), pad(BC0, 4)
            pa = procaware.iloc[0]
            fit = pad(float(pa.best_D_coefficient)*BD0 + float(pa.best_C_coefficient)*BC0, 4)
            fullnorm = fr.phase_RMS_full_mrad / phase_noise
            Fref = fit.copy()
            Fref[-1] = float(pa.two_process_relative_misfit) * fullnorm
        blocks[('SIP', t)] = block(BS, BD, BC, Fref, f'sip_{t}', t, 'SIP')
    return blocks


def assemble(blocks, modalities):
    Bs, Rs = [], []
    for mod in modalities:
        for t in [365, 1825]:
            b = blocks[(mod, t)]
            Bs.append(b['B'])
            Rs.append(b['R'])
    return np.vstack(Bs), np.concatenate(Rs)


def hpair(a):
    s, d, c = a[:, 0], a[:, 1], a[:, 2]
    return (s*d + s*c + d*c) / 3.0


def main():
    blocks = build_blocks(phase_noise=5.0)
    cases = {
        'seismic only': ['seismic'],
        'seismic + DC': ['seismic', 'DC'],
        'seismic + DC + SIP': ['seismic', 'DC', 'SIP'],
    }
    sob = qmc.Sobol(d=3, scramble=True, seed=1414)
    samples = 1.8 * sob.random_base2(18)
    h = hpair(samples)
    truth = np.array([[1., 1., 1.]])
    htruth = hpair(truth)[0]

    rows, weights = [], {}
    for name, mods in cases.items():
        B, R = assemble(blocks, mods)
        pred = samples @ B.T + h[:, None] * R[None, :]
        ytrue = (truth @ B.T).ravel() + htruth * R
        obs = ytrue + RNG.normal(0, 1, len(ytrue))
        w = normalize_logw(log_weights(pred, obs))
        weights[name] = w
        mean, std, cov, corr, q05, q50, q95 = weighted_stats(samples, w)
        rows.append([name, *mean, *std, corr[1, 2], q05[0], q95[0], q05[1], q95[1], q05[2], q95[2]])

    post = pd.DataFrame(rows, columns=['data_case','mean_aS','mean_aD','mean_aC','sd_aS','sd_aD','sd_aC','corr_D_C',
                                      'aS_p05','aS_p95','aD_p05','aD_p95','aC_p05','aC_p95'])
    post.to_csv(OUT / 'headline_posterior_reproduced.csv', index=False)

    # Null-direction variances
    null_rows = []
    uplus = np.array([1., 1.]) / np.sqrt(2)
    uminus = np.array([1., -1.]) / np.sqrt(2)
    for name in cases:
        _, _, cov, corr, *_ = weighted_stats(samples, weights[name])
        C = cov[1:3, 1:3]
        vp = float(uplus @ C @ uplus)
        vm = float(uminus @ C @ uminus)
        null_rows.append([name, np.sqrt(vp), np.sqrt(vm), corr[1, 2]])
    null = pd.DataFrame(null_rows, columns=['data_case','sd_DplusC','sd_DminusC','corr_D_C'])
    null.to_csv(OUT / 'null_direction_reproduced.csv', index=False)

    # Figure matching the manuscript's posterior comparison logic
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True, sharey=True)
    for ax, name in zip(axes, cases):
        w = weights[name]
        H, _, _ = np.histogram2d(samples[:, 1], samples[:, 2], bins=60,
                                 range=[[0, 1.8], [0, 1.8]], weights=w)
        ax.imshow(H.T, origin='lower', extent=[0,1.8,0,1.8], aspect='auto')
        ax.scatter([1], [1], marker='x', s=60)
        ax.set_title(name)
        ax.set_xlabel('aD')
    axes[0].set_ylabel('aC')
    fig.tight_layout()
    fig.savefig(OUT / 'posterior_DC_reproduced.png', dpi=190)
    plt.close(fig)

    # Verification against frozen headline tolerances
    frozen = pd.read_csv(DATA / 'basalt_co2_v1_4_headline_posterior.csv').set_index('data_case')
    repro = post.set_index('data_case')
    checks = []
    for case in frozen.index:
        for col in ['sd_aD','sd_aC','corr_D_C']:
            checks.append((case, col, float(repro.loc[case,col]), float(frozen.loc[case,col])))
    verify = pd.DataFrame(checks, columns=['data_case','metric','reproduced','frozen'])
    verify['abs_difference'] = np.abs(verify.reproduced - verify.frozen)
    verify.to_csv(OUT / 'verification.csv', index=False)
    print(post.to_string(index=False))
    print('\nmax absolute difference from frozen headline metrics:', verify.abs_difference.max())


if __name__ == '__main__':
    main()
