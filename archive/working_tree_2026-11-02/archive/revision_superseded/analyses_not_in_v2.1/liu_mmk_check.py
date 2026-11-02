"""Port of Liu et al. 03-Mann_Kendall_Trend_Test.m and 06-Modified-Mann-Kendall-Trend-Test.m;
false-positive rate on trend-free autocorrelated series (n = 372 monthly values)."""
import numpy as np
from scipy.stats import norm
import pymannkendall as mk

def mk_liu(x):
    n = len(x); S = sum(np.sign(x[k+1:] - x[k]).sum() for k in range(n-1))
    var = n*(n-1)*(2*n+5)/18
    Z = (S-1)/np.sqrt(var) if S > 0 else ((S+1)/np.sqrt(var) if S < 0 else 0)
    return 2*(1-norm.cdf(abs(Z))), S

def acf(x, L):  # MATLAB autocorr: biased, demeaned
    x = x - x.mean(); c0 = np.dot(x, x)/len(x)
    return np.array([np.dot(x[:len(x)-k], x[k:])/len(x)/c0 for k in range(L+1)])

def mmk_liu(x):
    n = len(x); _, S = mk_liu(x)
    L = min(24, n//4); r = acf(x, L)[1:]
    tau = 1 + 2*np.sum((1-np.arange(1, L+1)/n)*r)
    ne = n/tau; var = ne*(ne-1)*(2*ne+5)/18
    Z = (S-1)/np.sqrt(var) if S > 0 else ((S+1)/np.sqrt(var) if S < 0 else 0)
    return 2*(1-norm.cdf(abs(Z))), tau

rng = np.random.default_rng(1)
for phi, smooth in ((0.0, 1), (0.5, 1), (0.5, 4), (0.8, 4)):
    fp_mk = fp_mmk = fp_hr = 0; taus = []
    N = 300
    for _ in range(N):
        e = rng.standard_normal(372 + 50); x = np.zeros_like(e)
        for t in range(1, len(e)): x[t] = phi*x[t-1] + e[t]
        x = np.convolve(x, np.ones(smooth)/smooth, mode="valid")[-372:]
        p1, _ = mk_liu(x); p2, tau = mmk_liu(x); taus.append(tau)
        p3 = mk.hamed_rao_modification_test(x).p
        fp_mk += p1 < 0.05; fp_mmk += p2 < 0.05; fp_hr += p3 < 0.05
    print(f"AR(1) phi={phi}, {smooth}-month running mean: false-positive rate at 5%: "
          f"MK {fp_mk/N:.2f}, Liu 'modified MK' {fp_mmk/N:.2f}, pymannkendall Hamed-Rao {fp_hr/N:.2f}; median tau {np.median(taus):.1f}")
