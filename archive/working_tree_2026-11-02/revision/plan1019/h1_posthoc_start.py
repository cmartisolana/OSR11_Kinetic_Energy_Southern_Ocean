"""POST-HOC sensitivity (not part of the pre-specified gate): G1 coverage when the branch start is restricted to
climatological envelopes inside each site's ACC latitude band. Coverage only; no trends."""
import sys; sys.argv=[sys.argv[0], sys.argv[1], sys.argv[2], 'gate']
import numpy as np
exec(open('h1_halfpower_meanders.py').read().split('if __name__ == "__main__":')[0])
def branch_band(ccnt, ctab, cols, cov, band):
    best, start = -np.inf, None
    for q, j in enumerate(cols):
        for c in range(ccnt[j]):
            lh = ctab[j, c, FI["lat_hp"]]
            if band[0] <= lh <= band[1] and ctab[j, c, FI["peak_val"]] > best and cov[q, c] >= 0.9:
                best, start = ctab[j, c, FI["peak_val"]], (q, c)
    if start is None: return []
    path = {start[0]: start[1]}
    for step in (-1, 1):
        q, last = start[0], ctab[cols[start[0]], start[1], FI["lat_hp"]]; gap = 0
        while 0 <= q+step < len(cols):
            q += step; j = cols[q]
            cand = [(abs(ctab[j, c, FI["lat_hp"]]-last), c) for c in range(ccnt[j]) if abs(ctab[j, c, FI["lat_hp"]]-last) <= 0.25]
            if cand: c = min(cand)[1]; path[q] = c; last = ctab[j, c, FI["lat_hp"]]; gap = 0
            else:
                gap += 1
                if gap > 4: break
    return sorted(path.items())
for name in ["KERG", "PAR", "DRAKE"]:
    c = SITES[name]; lat, lon, U, V, S2, ICE = load_site(c); KE = 0.5*(U.astype(float)**2+V.astype(float)**2)
    cols = np.where((lon >= c["mea"][0]) & (lon <= c["mea"][1]))[0]
    ccnt, ctab = F.detect_field(np.nanmean(KE, 0), lat, F.meridional_dists(lat), THR, MINW, FRAC, MAXENV)
    cnt, tab = detect_stack(KE, lat); M = match(cnt, tab, ccnt, ctab, lat, cols); cov = np.isfinite(M).mean(0)
    br = branch_band(ccnt, ctab, cols, cov, c["band"])
    usable = sum(cov[q, cc] >= 0.9 for q, cc in br)
    anycov = np.mean([np.nanmax(np.where(np.isfinite(cov[q]), cov[q], 0)) >= 0.9 for q in range(len(cols))])
    print(name, "band-start branch usable fraction", round(usable/len(cols), 3), "| columns with ANY envelope cov>=0.9:", round(anycov, 3), flush=True)
