"""Replace container-specific absolute paths in the repository scripts by paths relative to revision/ (26 Oct 2026)."""
import re, sys
from pathlib import Path
ROOT = Path(sys.argv[1])  # revision dir
PREF = [("/mnt/user-data/outputs/revision/outputs", "_REV_OUT"), ("/mnt/user-data/uploads/repository/revision", "_REV"),
        ("/mnt/user-data/uploads/repository/outputs", "_ROOT_OUT"), ("/mnt/user-data/outputs/revision", "_REV"),
        ("/home/claude/work/drake/out", "_DK_OUT"), ("/home/claude/work/drake", "_DK"), ("/home/claude/w/cache", "_CACHE"),
        ("/home/claude/ext/TrendsDrakePassageTransp", "_GV")]
HDR = ("import os as _os\nfrom pathlib import Path as _P\n_REV = str(_P(__file__).resolve().parents[{k}])  # revision/\n"
       "_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'\n"
       "_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')\n")
lit = re.compile(r'''(?P<f>[fF]?)(?P<q>"|')(?P<s>[^"'\n]*?)(?P=q)''')
n_files = 0
for f in sorted(ROOT.rglob("*.py")):
    if "archive" in f.parts: continue
    s = f.read_text(encoding="utf8")
    if not any(p in s for p, _ in PREF): continue
    def sub(m):
        s0 = m.group("s")
        for p, v in PREF:
            if s0.startswith(p):
                rest = s0[len(p):]
                if m.group("f"): return f'({v} + f"{rest}")'
                return f'({v} + "{rest}")' if rest else v
        return m.group(0)
    new = lit.sub(sub, s)
    k = len(f.relative_to(ROOT).parts) - 1
    # insert header after the module docstring (if any) and __future__ imports
    m = re.match(r'(\s*(?:"""[\s\S]*?"""|\'\'\'[\s\S]*?\'\'\')\s*\n)?((?:from __future__[^\n]*\n)*)', new)
    pos = m.end() if m else 0
    new = new[:pos] + HDR.format(k=k) + new[pos:]
    f.write_text(new, encoding="utf8"); n_files += 1; print("patched", f.relative_to(ROOT), "depth", k)
print(n_files, "files")
