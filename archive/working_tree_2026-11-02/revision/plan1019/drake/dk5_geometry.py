import scipy.io as sio
L_ = lambda f, k: sio.loadmat(GV / 'Datasets' / f, squeeze_me=True, struct_as_record=False)[k]
R_E = 6.371e6
L = lambda f, k: sio.loadmat(GV / "Datasets" / f, squeeze_me=True, struct_as_record=False)[k]
TT = L_("TotalTransport_os38nbmaxz_760_2.mat", "totaltransp")
blat = np.array(TT.lat, float); blon = np.array(TT.lon, float) - 360; dist = np.array(TT.dist, float)
keep = dist <= 840
mlat = np.nanmean(blat[keep], 1); mlon = np.nanmean(blon[keep], 1)
s = np.r_[0, np.cumsum(np.hypot(np.diff(mlon) * np.cos(np.radians(mlat[:-1])), np.diff(mlat)))]
ss = np.linspace(0, s[-1], 121); lo = np.interp(ss, s, mlon); la = np.interp(ss, s, mlat)
dkm = np.interp(ss, s, dist[keep])                              # distance from north (km), GV23 convention
latm = np.radians(0.5 * (la[:-1] + la[1:]))
seg = np.hypot(R_E * np.cos(latm) * np.radians(np.diff(lo)), R_E * np.radians(np.diff(la)))
dsx = np.r_[seg[0] / 2, 0.5 * (seg[:-1] + seg[1:]), seg[-1] / 2]
tx = np.gradient(lo) * np.cos(np.radians(la)); ty = np.gradient(la); nn = np.hypot(tx, ty)
nx, ny = -ty / nn, tx / nn; flip = nx < 0; nx[flip] *= -1; ny[flip] *= -1
