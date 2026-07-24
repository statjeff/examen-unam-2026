#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modelo jerarquico bayesiano de cuantil censurado para los puntajes de corte
del concurso de seleccion a licenciatura de la UNAM, 2021-2026.

Estructura:
  - Los aciertos minimos son el estadistico de orden K-esimo de la distribucion
    latente de puntajes. Se invierten hacia la media latente mu_jt.
  - Modelo jerarquico: mu_jt = a_j + tau_t + eta_jt, con a_j efecto de
    carrera-plantel-modalidad y tau_t caminata aleatoria de primer orden.
  - Los ~200 efectos a_j se MARGINALIZAN ANALITICAMENTE (identidad de
    Sherman-Morrison), dejando ~12 hiperparametros globales para emcee.
  - Validacion leave-one-year-out (LOYO) para calibrar la discrepancia
    predictiva esperada bajo el modelo sin choque.

Autor del analisis: Jeffrey E. Barcenas Mosqueda
"""

import numpy as np
import pandas as pd
from scipy.stats import norm
import emcee
import os, json

rng = np.random.default_rng(20260721)
OUT = "/home/claude/salidas"
os.makedirs(OUT, exist_ok=True)

PISO = 40          # piso administrativo de aciertos (censura por la izquierda)
N_ITEMS = 120      # reactivos del examen

# ----------------------------------------------------------------------
# 1. Preparacion de datos
# ----------------------------------------------------------------------

def cargar():
    d = pd.read_csv("/mnt/user-data/uploads/metadata_carreras.csv")
    # La llave del panel es (codigo, modalidad): un mismo codigo se reutiliza
    # para las modalidades 'abierta' y 'suayed'.
    d["unidad"] = d.codigo.astype(str) + "_" + d.modalidad
    d["p"] = d.seleccionados / d.presentaron_examen

    # Criterio de inclusion. Se descartan observaciones donde el corte no es
    # informativo sobre la distribucion latente:
    #   (i)  corte en el piso administrativo -> censurado por la izquierda
    #   (ii) saturacion: practicamente todo sustentante fue admitido
    #   (iii) cola extrema: p muy chico hace explotar la varianza de mu_hat
    d["usable"] = (
        (d.aciertos_minimos > PISO)
        & (d.p < 0.95)
        & (d.p > 0.01)
        & (d.presentaron_examen > 30)
    )
    return d


def momentos_mu(d, s):
    """Inversion del cuantil: dada la desviacion estandar latente s, devuelve
    la estimacion puntual de la media latente y su varianza por delta method.

        p_jt = P(X >= c_jt),  z_jt = Phi^{-1}(1 - p_jt),  mu_jt = c_jt - s*z_jt

    La incertidumbre proviene del error binomial en p_hat = K/N.
    """
    p = d["p"].values
    N = d["presentaron_examen"].values
    c = d["aciertos_minimos"].values.astype(float)
    z = norm.ppf(1.0 - p)
    mu_hat = c - 0.5 - s * z          # -0.5: correccion por continuidad
    var_p = p * (1.0 - p) / N
    dens = norm.pdf(z)
    var_mu = (s ** 2) * var_p / np.maximum(dens ** 2, 1e-12)
    return mu_hat, var_mu, z


# ----------------------------------------------------------------------
# 2. Verosimilitud marginalizada
# ----------------------------------------------------------------------

class Modelo:
    """Verosimilitud con los efectos de unidad a_j integrados analiticamente.

    Para la unidad j con observaciones t en su ventana:
        r_jt = mu_hat_jt - tau_t
        r_j ~ MVN( m_{area(j)} * 1 ,  D_j + sigma_a^2 * J )
    con D_j = diag(v_jt + sigma_eta^2). Sherman-Morrison da la inversa y el
    determinante en O(n_j).
    """

    def __init__(self, d, anios):
        self.anios = list(anios)
        sub = d[d.year.isin(self.anios) & d.usable].copy()
        # conservar unidades con al menos 3 anios utilizables en la ventana
        cnt = sub.groupby("unidad").size()
        minobs = min(3, len(self.anios))
        keep = cnt[cnt >= minobs].index
        sub = sub[sub.unidad.isin(keep)].copy()

        self.d = sub
        self.unidades = sorted(sub.unidad.unique())
        self.uidx = {u: i for i, u in enumerate(self.unidades)}
        self.J = len(self.unidades)
        self.tidx = {y: k for k, y in enumerate(self.anios)}

        self.ju = sub.unidad.map(self.uidx).values.astype(int)
        self.jt = sub.year.map(self.tidx).values.astype(int)
        self.area = sub.area.values.astype(int) - 1
        # area por unidad
        self.area_u = np.zeros(self.J, dtype=int)
        for u, a in zip(self.ju, self.area):
            self.area_u[u] = a
        self.c = sub.aciertos_minimos.values.astype(float)
        self.p = sub.p.values
        self.N = sub.presentaron_examen.values.astype(float)
        self.T = len(self.anios)

        # --- precomputo para la inversion del cuantil (independiente de s) ---
        self.z = norm.ppf(1.0 - self.p)
        self.zbar = float(self.z.mean())   # centrado: decorrela s de m
        self.zc = self.z - self.zbar
        self.dens2 = np.maximum(norm.pdf(self.z) ** 2, 1e-12)
        self.varp = self.p * (1.0 - self.p) / self.N
        self.ratio = self.varp / self.dens2      # var_mu = s^2 * ratio
        self.area_obs = self.area_u[self.ju]
        self.n_j = np.bincount(self.ju, minlength=self.J).astype(float)

    # --- empaquetado de parametros ---
    # theta = [tau_2..tau_T (T-1), m_1..m_4 (4), log s, log sig_a, log sig_eta, log sig_rw]
    @property
    def ndim(self):
        return (self.T - 1) + 4 + 4

    def desempaca(self, th):
        T = self.T
        tau = np.concatenate([[0.0], th[: T - 1]])
        m = th[T - 1 : T + 3]
        s, sa, se, srw = np.exp(th[T + 3 : T + 7])
        return tau, m, s, sa, se, srw

    def log_prior(self, th):
        tau, m, s, sa, se, srw = self.desempaca(th)
        lp = 0.0
        # s: desviacion estandar latente de puntajes (aciertos sobre 120).
        # Prior informativo: es el parametro debilmente identificado.
        lp += norm.logpdf(np.log(s), np.log(17.0), 0.25)
        # semi-normales sobre las escalas jerarquicas (+ jacobiano log)
        for v, sc in ((sa, 20.0), (se, 6.0), (srw, 6.0)):
            if v <= 0 or v > 200:
                return -np.inf
            lp += norm.logpdf(v, 0.0, sc) + np.log(v)
        lp += np.sum(norm.logpdf(m, 60.0, 30.0))
        # caminata aleatoria de primer orden sobre los efectos de anio
        lp += np.sum(norm.logpdf(np.diff(tau), 0.0, srw))
        if not np.isfinite(lp):
            return -np.inf
        return lp

    def log_like(self, th):
        """Verosimilitud con a_j marginalizado, vectorizada sobre unidades.

        Sherman-Morrison aplicado en bloque: para Sigma_j = D_j + sa^2 J,
            quad_j   = S(r^2 u) - sa^2 S(ru)^2 / (1 + sa^2 S(u))
            logdet_j = -S(log u) + log1p(sa^2 S(u))
        donde cada S(.) es una suma por unidad calculada con bincount.
        """
        tau, m, s, sa, se, srw = self.desempaca(th)
        mu_hat = self.c - 0.5 - s * self.zc
        var_mu = (s ** 2) * self.ratio
        r = mu_hat - tau[self.jt] - m[self.area_obs]
        u = 1.0 / (var_mu + se ** 2)

        J, ju = self.J, self.ju
        Su = np.bincount(ju, weights=u, minlength=J)
        Sru = np.bincount(ju, weights=r * u, minlength=J)
        Sr2u = np.bincount(ju, weights=r * r * u, minlength=J)
        Slu = np.bincount(ju, weights=np.log(u), minlength=J)

        den = 1.0 + sa ** 2 * Su
        quad = Sr2u - (sa ** 2) * Sru ** 2 / den
        logdet = -Slu + np.log(den)
        return -0.5 * float(np.sum(quad + logdet + self.n_j * np.log(2 * np.pi)))

    def log_post(self, th):
        lp = self.log_prior(th)
        if not np.isfinite(lp):
            return -np.inf
        ll = self.log_like(th)
        return lp + ll if np.isfinite(ll) else -np.inf

    def inicio(self, nw):
        """Arranque guiado por los datos. emcee con movimientos diferenciales
        necesita que el ensamble inicial tenga la escala del posterior en cada
        direccion; una bola isotropa minuscula mezcla pesimo."""
        s0 = 17.0
        mu0 = self.c - 0.5 - s0 * self.zc
        # niveles por area
        m0 = np.array([mu0[self.area_obs == g].mean() if np.any(self.area_obs == g)
                       else mu0.mean() for g in range(4)])
        # escalas jerarquicas empiricas
        res = mu0 - m0[self.area_obs]
        med_u = np.bincount(self.ju, weights=res, minlength=self.J) / self.n_j
        sa0 = max(med_u.std(), 1.0)
        se0 = max((res - med_u[self.ju]).std(), 0.5)
        th0 = np.concatenate([np.zeros(self.T - 1), m0,
                              np.log([s0, sa0, se0, 1.5])])
        # dispersion por bloque: tau ~1, m ~2, log-escalas ~0.15
        esc = np.concatenate([np.full(self.T - 1, 1.0), np.full(4, 2.0),
                              np.full(4, 0.15)])
        return th0 + esc * rng.normal(size=(nw, self.ndim))


def ajustar(modelo, nsteps=4000, nwalkers=None, etiqueta=""):
    nd = modelo.ndim
    nw = nwalkers or max(4 * nd, 40)
    moves = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)]
    sam = emcee.EnsembleSampler(nw, nd, modelo.log_post, moves=moves)
    sam.run_mcmc(modelo.inicio(nw), nsteps, progress=False)
    try:
        tau_int = sam.get_autocorr_time(quiet=True)
        burn = int(3 * np.nanmax(tau_int))
        thin = max(1, int(np.nanmax(tau_int) / 2))
    except Exception:
        burn, thin = nsteps // 3, 10
    burn = min(max(burn, nsteps // 4), nsteps // 2)
    ch = sam.get_chain(discard=burn, thin=thin, flat=True)
    af = np.mean(sam.acceptance_fraction)
    print(f"  [{etiqueta}] ndim={nd} nw={nw} burn={burn} thin={thin} "
          f"muestras={len(ch)} aceptacion={af:.3f} tau_max={np.nanmax(tau_int):.1f}")
    return ch, sam, {"aceptacion": float(af), "tau_max": float(np.nanmax(tau_int)),
                     "n_muestras": int(len(ch)), "burn": int(burn)}


# ----------------------------------------------------------------------
# 3. Validacion leave-one-year-out
# ----------------------------------------------------------------------

def loyo(d, anio_test, nsteps=7000, ndraw=800):
    """Entrena con los anios previos, predice el corte del anio excluido
    manteniendo fija la tasa de seleccion observada K/N."""
    anios_tr = [y for y in sorted(d.year.unique()) if y < anio_test]
    if len(anios_tr) < 2:
        return None
    mod = Modelo(d, anios_tr)
    ch, _, diag = ajustar(mod, nsteps=nsteps, nwalkers=80, etiqueta=f"LOYO {anio_test}")

    # unidades evaluables: usables en el anio de prueba y presentes en el ajuste
    te = d[(d.year == anio_test) & d.usable].copy()
    te = te[te.unidad.isin(mod.unidades)].copy()

    idx = rng.choice(len(ch), size=min(ndraw, len(ch)), replace=False)
    D = np.zeros((len(idx), len(te)))          # discrepancia c_obs - c_pred
    Cpred = np.zeros_like(D)

    for k, ii in enumerate(idx):
        tau, m, s, sa, se, srw = mod.desempaca(ch[ii])
        # extrapolacion de la caminata aleatoria al anio de prueba
        tau_new = tau[-1] + rng.normal(0.0, srw)

        # posterior conjugada de a_j dados los datos de entrenamiento
        mu_hat = mod.c - 0.5 - s * mod.zc
        var_mu = (s ** 2) * mod.ratio
        r = mu_hat - tau[mod.jt] - m[mod.area_obs]
        w = 1.0 / (var_mu + se ** 2)
        Su = np.bincount(mod.ju, weights=w, minlength=mod.J)
        Sru = np.bincount(mod.ju, weights=w * r, minlength=mod.J)
        prec = 1.0 / sa ** 2 + Su
        a_mean = Sru / prec
        a_draw = a_mean + rng.normal(size=mod.J) / np.sqrt(prec)

        ju_te = te.unidad.map(mod.uidx).values
        area_te = te.area.values.astype(int) - 1
        mu_pred = (a_draw[ju_te] + m[area_te] + tau_new
                   + rng.normal(0.0, se, size=len(te)))
        z_te = norm.ppf(1.0 - te.p.values) - mod.zbar   # mismo centrado
        c_pred = mu_pred + s * z_te + 0.5
        Cpred[k] = c_pred
        D[k] = te.aciertos_minimos.values - c_pred

    return {"anio": anio_test, "te": te, "D": D, "Cpred": Cpred,
            "n_unidades": len(te), "diag": diag, "modelo": mod, "chain": ch}


# ----------------------------------------------------------------------
# 4. Ejecucion
# ----------------------------------------------------------------------

if __name__ == "__main__":
    d = cargar()
    print(f"Observaciones totales: {len(d)} | usables: {d.usable.sum()} "
          f"({100*d.usable.mean():.1f}%)")
    print(d.groupby('year').usable.agg(['sum', 'size']))

    # --- 4a. Ajuste completo 2021-2026 (tau_2026 libre) ---
    print("\n=== Ajuste completo 2021-2026 ===")
    mod_full = Modelo(d, sorted(d.year.unique()))
    ch_full, sam_full, diag_full = ajustar(mod_full, nsteps=8000, nwalkers=80, etiqueta="completo")
    np.save(f"{OUT}/chain_full.npy", ch_full)
    with open(f"{OUT}/unidades_full.json", "w") as f:
        json.dump({"unidades": list(map(str, mod_full.unidades)),
                   "anios": list(map(int, mod_full.anios)),
                   "diag": diag_full, "zbar": mod_full.zbar}, f)

    tau_post = np.array([mod_full.desempaca(t)[0] for t in ch_full])
    s_post = np.array([mod_full.desempaca(t)[2] for t in ch_full])
    print("\nEfecto latente de anio tau_t (aciertos, base 2021=0):")
    for k, y in enumerate(mod_full.anios):
        q = np.percentile(tau_post[:, k], [2.5, 16, 50, 84, 97.5])
        print(f"  {y}: mediana {q[2]:7.2f}  IC68 [{q[1]:6.2f},{q[3]:6.2f}]  "
              f"IC95 [{q[0]:6.2f},{q[4]:6.2f}]")
    dif = tau_post[:, -1] - tau_post[:, -2]
    print(f"\n  Delta tau(2026-2025): mediana {np.median(dif):.2f}  "
          f"IC95 [{np.percentile(dif,2.5):.2f},{np.percentile(dif,97.5):.2f}]  "
          f"P(delta>0)={np.mean(dif>0):.4f}")
    print(f"  s (DE latente de puntajes): {np.median(s_post):.2f} "
          f"[{np.percentile(s_post,2.5):.2f},{np.percentile(s_post,97.5):.2f}]")

    # --- 4b. Validacion LOYO ---
    print("\n=== Validacion leave-one-year-out ===")
    res = {}
    for y in [2023, 2024, 2025, 2026]:
        r = loyo(d, y, nsteps=7000, ndraw=800)
        if r is None:
            continue
        Dm = r["D"].mean(axis=1)          # discrepancia media por draw
        q = np.percentile(Dm, [2.5, 50, 97.5])
        rmse = np.sqrt((r["D"] ** 2).mean())
        print(f"  {y}: n={r['n_unidades']:3d}  D_media mediana={q[1]:6.2f} "
              f"IC95 [{q[0]:6.2f},{q[2]:6.2f}]  RMSE={rmse:5.2f}")
        res[y] = {"Dm": Dm, "D": r["D"], "te": r["te"], "n": r["n_unidades"],
                  "rmse": float(rmse), "diag": r["diag"]}
        np.save(f"{OUT}/D_{y}.npy", r["D"])
        r["te"].to_csv(f"{OUT}/te_{y}.csv", index=False)

    # --- 4c. Heterogeneidad de la discrepancia en 2026 ---
    if 2026 in res:
        te26 = res[2026]["te"].copy()
        te26["D_med"] = np.median(res[2026]["D"], axis=0)
        print("\n=== Discrepancia 2026 por area ===")
        print(te26.groupby("area").D_med.agg(["mean", "median", "count"]).round(2))
        print("\n=== Discrepancia 2026 por modalidad ===")
        print(te26.groupby("modalidad").D_med.agg(["mean", "median", "count"]).round(2))
        te26.to_csv(f"{OUT}/discrepancia_2026.csv", index=False)

    # --- 4d. Sensibilidad al prior sobre s ---
    print("\n=== Sensibilidad: s fijo ===")
    sens = {}
    for s_fix in [13.0, 17.0, 21.0]:
        class ModFix(Modelo):
            _s = s_fix
            def desempaca(self, th):
                tau, m, s, sa, se, srw = Modelo.desempaca(self, th)
                return tau, m, self._s, sa, se, srw
        mf = ModFix(d, sorted(d.year.unique()))
        chf, _, _ = ajustar(mf, nsteps=5000, nwalkers=80, etiqueta=f"s={s_fix}")
        tp = np.array([mf.desempaca(t)[0] for t in chf])
        dd = tp[:, -1] - tp[:, -2]
        sens[s_fix] = np.percentile(dd, [2.5, 50, 97.5])
        print(f"  s={s_fix}: Delta tau(2026-2025) = {sens[s_fix][1]:.2f} "
              f"[{sens[s_fix][0]:.2f},{sens[s_fix][2]:.2f}]")

    np.save(f"{OUT}/tau_post_full.npy", tau_post)
    np.save(f"{OUT}/s_post_full.npy", s_post)
    with open(f"{OUT}/sensibilidad.json", "w") as f:
        json.dump({str(k): list(map(float, v)) for k, v in sens.items()}, f)
    print("\nListo. Salidas en", OUT)
