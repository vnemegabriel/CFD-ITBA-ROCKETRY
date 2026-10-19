"""python -m finflutter <command> ... ; python -m finflutter <case.yaml> is the same as 'run'."""
import argparse
import json
import re
from pathlib import Path

import numpy as np

from . import aero, aeroelastic, cfd, flight, gaf, lco, plots, static, taptest
from .run import build_structure, load_case, main as run_main


def semichord(fin):
    return (fin.cr + fin.ct) / 4


def loss_factors(c, n, g_file=None):
    if g_file:
        d = json.loads(Path(g_file).read_text())
        g = np.array([m.get("g_mean", np.nan) for m in d["modes"]], float)
        g = np.concatenate([g, np.full(max(0, n - len(g)), np.nan)])[:n]
        return np.where(np.isfinite(g), g, np.nanmin(g) if np.isfinite(g).any() else c["g"])
    return np.broadcast_to(np.asarray(c["g"], float), (n,)).copy()


def cmd_run(a):
    run_main(a.case)


def cmd_figures(a):
    s = build_structure(load_case(a.case))
    out = Path(a.out)
    print(plots.geometry(s["fin"], out / "geometry.png"))
    print(plots.modes(s["nodes"], s["elems"], s["phi"], s["omega"], out / "modes.png"))
    c = load_case(a.case)
    rows = []
    for nx, ny in ((10, 4), (20, 8), (40, 16), (80, 32)):
        c["mesh"], c["modes"] = {"nx": nx, "ny": ny}, 4
        rows.append((f"{nx}x{ny}", nx * ny, list(build_structure(c)["omega"] / 2 / np.pi)))
    print(plots.convergence(rows, out / "fe_convergence.png"))


def cmd_modes(a):
    s = build_structure(load_case(a.case))
    n, maxabs = s["phi"].shape[1], np.abs(s["phi"]).max(axis=0)
    if a.signal == "pulse":
        sig = cfd.pulse(n, a.mode, a.amplitude, maxabs, a.fmax)
    elif a.signal == "harmonic":
        crit = json.loads(Path(a.critical).read_text())
        v = np.array(crit["vector_re"]) + 1j * np.array(crit["vector_im"])
        sig = cfd.harmonic(s["phi"], v, a.amplitude, 2 * np.pi * crit["f_Hz"])
    else:
        sig = cfd.no_signal(n)
    roll = 45.0 if a.sector == "half" else 0.0
    if a.run:
        sector = cfd.run_sector(a.run)
        roll = 45.0 if sector == "half" else 0.0
        if sector != "full":
            print(f"!! {a.run} is a {sector} mesh: the mirror image of fin 0 moves with it")
        path = Path(a.run) / "constant" / "flutterModes"
    else:
        path = Path(a.out)
    cfd.write_flutter_modes(path, s, sig, roll_deg=roll)
    print(f"{path}: {n} modes, f = {' '.join(f'{w/2/np.pi:.0f}' for w in s['omega'])} Hz, signal {sig['type']}")
    if sig["type"] == "pulse":
        print(f"  mode {a.mode}, tip amplitude {a.amplitude*1e3:.2f} mm, eta peak {sig['coef'][a.mode-1]:.4g}, "
              f"tau {sig['tau']*1e3:.3f} ms, t0 {sig['t0']*1e3:.3f} ms")


def cmd_conditions(a):
    c = load_case(a.case)
    pts = flight.load(Path(c["flight"]["csv"]), 0.0)
    asc = np.arange(np.argmax(pts["mach"]) + 1)
    print(f"{'Mach':>6} {'t [s]':>7} {'h [m]':>7} {'p [Pa]':>8} {'T [K]':>7} {'U [m/s]':>8} {'q [Pa]':>8}")
    for M in a.mach:
        i = asc[np.argmin(np.abs(pts["mach"][asc] - M))]
        rho = pts["rho"][i]
        T = (pts["U"][i] / pts["mach"][i]) ** 2 / (1.4 * 287.053)
        p = rho * 287.053 * T
        print(f"{pts['mach'][i]:6.3f} {pts['t'][i]:7.2f} {pts['h'][i]:7.0f} {p:8.0f} {T:7.2f} {pts['U'][i]:8.1f} {pts['q'][i]:8.0f}"
              f"   --Minf {M:g} --pInf {p:.0f} --Tinf {T:.2f}")


def cmd_static(a):
    s = build_structure(load_case(a.case))
    r = static.compare(a.run_alpha, a.run_zero, s)
    print(f"M {r['mach']:.2f}, alpha {r['alpha_deg']:.1f} deg -> incidence of fin 0 {r['theta_deg']:.2f} deg, q {r['q']:.0f} Pa")
    print(f"{'mode':>4} {'F_cfd':>10} {'F_piston':>10} {'cfd/piston':>10}")
    for k in range(min(a.n, len(r["F_cfd"]))):
        print(f"{k+1:4d} {r['F_cfd'][k]:10.4g} {r['F_piston'][k]:10.4g} {r['F_cfd'][k]/r['F_piston'][k]:10.3f}")


def cmd_gaf(a):
    s = build_structure(load_case(a.case))
    d = gaf.assemble(a.runs, a.zero, b=semichord(s["fin"]))
    np.savez(a.out, **d)
    print(f"{a.out}: M {float(d['mach']):.3f}, k 0..{d['k'][-1]:.3f}, modes {list(d['modes'] + 1)}"
          + (f", missing {list(d['missing'] + 1)}" if len(d["missing"]) else ""))
    if a.plot:
        print(plots.gaf(d, a.plot))


def cmd_linearity(a):
    j, f, err = gaf.linearity(a.run_a, a.run_b, a.zero)
    band = f <= a.fmax
    print(f"mode {j+1}: max |Q(2A) - Q(A)| / max|Q| = {err[band].max():.3%} for f <= {a.fmax:g} Hz")
    if a.plot:
        print(plots.linearity(f, err, j, a.plot))


def cmd_flutter(a):
    c = load_case(a.case)
    s = build_structure(c)
    omega = s["omega"][: a.modes]
    g = loss_factors(c, len(omega), a.g_file)
    table = aeroelastic.RogerTable(gaf.table(a.gafs))
    for M, f in zip(table.machs, table.fits):
        print(f"M {M:.3f}: Roger fit residual {f['residual']:.2%}")
    pts = flight.load(Path(c["flight"]["csv"]), c["flight"].get("mach_min", 1.05))
    env = aeroelastic.envelope_rfa(omega, g, table, pts, c.get("q_factor", 20.0))
    if not env:
        raise SystemExit("no flight point inside the Mach range of the GAF files")
    crit = min(env, key=lambda r: r["margin"])
    rfa = table(crit["mach"])
    lam, v = aeroelastic.least_stable(omega, g, rfa, crit["q"], crit["U"])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "flutter.json").write_text(json.dumps({"g": g.tolist(), "envelope": env}, indent=1, default=float))
    critical = {"mach": crit["mach"], "q": crit["q"], "U": crit["U"], "margin": crit["margin"], "kind": crit["kind"],
                "f_Hz": float(abs(lam.imag) / 2 / np.pi), "sigma": float(lam.real),
                "vector_re": v.real.tolist(), "vector_im": v.imag.tolist()}
    (out / "critical.json").write_text(json.dumps(critical, indent=1, default=float))
    print(f"min margin q_crit/q = {crit['margin']:.2f} ({crit['kind'] or 'stable'}) at M = {crit['mach']:.3f}, "
          f"q = {crit['q']:.0f} Pa; least stable root there: {critical['f_Hz']:.0f} Hz, "
          f"damping ratio {-lam.real/abs(lam):.4f}")
    pt = aero.PistonTheory(s["phi"][:, : a.modes], s["mats"], theory="ackeret")
    bound = aeroelastic.envelope(omega, g, pt, pts, c.get("q_factor", 20.0))
    print(plots.envelope(env, bound, out / "margin.png"))
    print("->", out / "flutter.json", out / "critical.json")


def cmd_tap(a):
    s = build_structure(load_case(a.case))
    f_model = s["omega"][: a.modes] / 2 / np.pi
    rows = taptest.analyse(a.recording, n_modes=a.modes)
    summ = taptest.summary(rows, f_model)
    print(f"{'mode':>4} {'f model':>8} {'taps':>4} {'f meas':>8} {'+-':>5} {'g':>7} {'+-':>6} {'dg/damp':>8}")
    for m in summ:
        if m["n"]:
            print(f"{m['mode']:4d} {m['f_model']:8.0f} {m['n']:4d} {m['f_mean']:8.1f} {m['f_std']:5.1f} "
                  f"{m['g_mean']:7.4f} {m['g_std']:6.4f} {m['dg_damp']:8.4f}")
        else:
            print(f"{m['mode']:4d} {m['f_model']:8.0f}    0   (not found)")
    Path(a.out).write_text(json.dumps({"recording": str(a.recording), "taps": rows, "modes": summ}, indent=1, default=float))
    print("->", a.out)
    if a.plot:
        print(plots.tap(rows, f_model, a.plot))


def cmd_lco(a):
    c = load_case(a.case)
    s = build_structure(c)
    g = loss_factors(c, len(s["omega"]), a.g_file)
    rows, cross = lco.balance(a.runs, g, a.q_scale)
    print(f"{'A [mm]':>7} {'W_a [J]':>11} {'W_d [J]':>11} {'W_a/W_d':>8}")
    for r in rows:
        print(f"{r['A']*1e3:7.3f} {r['W_a']:11.4g} {r['W_d']:11.4g} {r['ratio']:8.3f}")
    for x in cross:
        print(f"LCO at A = {x['A']*1e3:.3f} mm ({'stable' if x['stable'] else 'unstable'})")
    if not cross:
        print("no crossing in the tested amplitudes")
    if a.plot:
        print(plots.lco(rows, cross, a.plot))


def cmd_sweep(a):
    """Forced-run sweep file and its flutterModes files, from the alpha-0 base runs of base sweep files."""
    c = load_case(a.case)
    s = build_structure(c)
    n, maxabs = s["phi"].shape[1], np.abs(s["phi"]).max(axis=0)
    mdir = Path(a.modes_dir).expanduser()
    mdir.mkdir(parents=True, exist_ok=True)

    def write(name, sig):
        cfd.write_flutter_modes(mdir / name, s, sig, roll_deg=45.0)
        return f"{a.modes_dir.rstrip('/')}/{name}"

    files = {"zero": write("none", cfd.no_signal(n))}
    for k in a.modes:
        for x in (1, 2):
            A = x * a.amplitude
            files[(k, x)] = write(f"pulse-m{k}-A{A*1e3:g}mm", cfd.pulse(n, k, A, maxabs, a.fmax))
    lines = ["# One run per line:  <name> <newCase.sh options>.  Written by",
             f"#     python -m finflutter sweep {' '.join(a.argv)}",
             "# Run it with ./runCase.sh [--np <N>] [--socket <S>] <this file>.",
             "# Forced-motion runs of the fin flutter analysis (fin_flutter/docs/PROCEDURE.md, step 7):",
             "# per Mach, A = 0 (noise of the base flow) and one pulse per mode; at the",
             f"# linearity Mach numbers also 2A for modes {a.modes[:2]}.  Each starts from its base run.", ""]
    for base in a.base:
        for line in Path(base).read_text().splitlines():
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            o = dict(zip(parts[1::1], parts[2::1]))
            if o.get("--solver") != a.base_solver or o.get("--alpha") != "0" or o.get("--mesh") != a.mesh:
                continue
            M = float(o["--Minf"])
            cond = f"--Minf {o['--Minf']} --pInf {o['--pInf']} --Tinf {o['--Tinf']}"
            tag = parts[0].replace("-pimple", "").replace("-a00-", "-").replace("m", "f", 1)
            common = f"--mesh {a.mesh} --solver forced --alpha 0  {cond}  --init {parts[0]}"
            lines.append(f"{tag+'-zero':<18} {common} --flutterModes {files['zero']}")
            for k in a.modes:
                lines.append(f"{tag+f'-p{k}':<18} {common} --flutterModes {files[(k, 1)]}")
            if any(abs(M - m) < 1e-6 for m in a.linearity_machs):
                for k in a.modes[:2]:
                    lines.append(f"{tag+f'-p{k}x2':<18} {common} --flutterModes {files[(k, 2)]}")
            lines.append("")
    Path(a.out).write_text("\n".join(lines))
    print(f"{a.out}: {sum(1 for l in lines if l and not l.startswith('#'))} runs; flutterModes files in {mdir}")


def check_forced(run, case=None):
    """Problems of a forced-motion run directory, before it runs. Returns a list of strings (empty = ready)."""
    run = Path(run)
    bad = []
    fm = run / "constant" / "flutterModes"
    if not fm.is_file():
        return [f"no {fm}: make the case with newCase.sh --flutterModes <file>"]
    text = fm.read_text()
    get = lambda k: re.search(rf"^\s*{k}\s+([^;]+);", text, re.M)
    roll = float(get("rollDeg").group(1))
    sector = cfd.run_sector(run)
    if (sector == "half") != (roll == 45.0):
        bad.append(f"rollDeg {roll:g} does not fit a {sector} mesh (half: 45, full: 0): write the file with --sector")
    ctrl = (run / "system" / "controlDict").read_text()
    app = re.search(r"^application\s+(\S+);", ctrl, re.M).group(1)
    if app != "rhoPimpleFoam":
        bad.append(f"application {app}: forced runs use rhoPimpleFoam (explicit solvers need dt ~1e-9 s on h0)")
    dt = float(re.search(r"^deltaT\s+([^;]+);", ctrl, re.M).group(1))
    end = float(re.search(r"^endTime\s+([^;]+);", ctrl, re.M).group(1))
    sig = text[text.index("signal"):]
    typ = re.search(r"type\s+(\w+);", sig).group(1)
    if typ == "pulse":
        tau = float(re.search(r"tau\s+([^;]+);", sig).group(1)); t0 = float(re.search(r"t0\s+([^;]+);", sig).group(1))
        if tau / dt < 50:
            bad.append(f"pulse tau {tau:g} s is only {tau/dt:.0f} time steps: the FFT needs >= 50 (smaller deltaT or larger tau)")
        if end < t0 + 4 * tau + 2e-3:
            bad.append(f"endTime {end:g} s ends before the response decays (t0 + 4 tau + 2 ms = {t0 + 4*tau + 2e-3:g} s)")
    if typ == "harmonic":
        w = float(re.search(r"omega\s+([^;]+);", sig).group(1)); tr = float(re.search(r"tRamp\s+([^;]+);", sig).group(1))
        if end < tr + 3 * 2 * np.pi / w:
            bad.append(f"endTime {end:g} s: the energy needs the ramp ({tr:g} s) plus 3 cycles")
    amp = re.search(r"amplitude\s+([^;]+);", sig)
    if amp and float(amp.group(1)) > 2e-3:
        bad.append(f"amplitude {float(amp.group(1))*1e3:g} mm: above 2 mm the mesh motion is not tested")
    init = (run / "initFrom")
    if not init.is_file():
        bad.append("no initFrom: a forced run starts from its base run (newCase.sh --init)")
    else:
        src = Path(init.read_text().strip())
        h = lambda f: re.search(r'meshHash\s+"(\w+)"', f.read_text()).group(1) if f.is_file() else None
        if h(src / "run.info") != h(run / "run.info"):
            bad.append(f"base {src} is on another mesh (meshHash differs): mapFields interpolates or crashes")
        if not any(p.name.replace(".", "").replace("e-", "").isdigit() and p.name != "0" for p in src.iterdir()):
            bad.append(f"base {src} has no solution time directory")
    return bad


def cmd_check(a):
    bad = check_forced(a.run)
    for b in bad:
        print("!!", b)
    print("ready" if not bad else f"{len(bad)} problem(s)")
    raise SystemExit(1 if bad else 0)


def main(argv):
    if argv and argv[0].endswith((".yaml", ".yml")):
        argv = ["run"] + argv
    p = argparse.ArgumentParser(prog="python -m finflutter")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="FE modes + piston-theory envelope (the bound)")
    r.add_argument("case"); r.set_defaults(f=cmd_run)

    r = sub.add_parser("figures", help="geometry, modes and FE convergence figures")
    r.add_argument("case"); r.add_argument("-o", "--out", default="docs/figures"); r.set_defaults(f=cmd_figures)

    r = sub.add_parser("modes", help="write constant/flutterModes for an OpenFOAM run")
    r.add_argument("case"); r.add_argument("run", nargs="?", help="run directory (writes <run>/constant/flutterModes)")
    r.add_argument("-o", "--out", default="flutterModes", help="file, when no run directory is given")
    r.add_argument("--signal", choices=["none", "pulse", "harmonic"], default="none")
    r.add_argument("--mode", type=int, default=1, help="pulse: the mode that moves (from 1)")
    r.add_argument("--amplitude", type=float, default=5e-4, help="peak displacement of the fin [m]")
    r.add_argument("--fmax", type=float, default=1500.0, help="pulse: frequency where its spectrum is 1 %% of the peak")
    r.add_argument("--critical", help="harmonic: critical.json written by 'flutter'")
    r.add_argument("--sector", choices=["half", "full"], default="half",
                   help="mesh the file is for, when no run directory is given: half = fins at 45 deg")
    r.set_defaults(f=cmd_modes)

    r = sub.add_parser("conditions", help="flight state at given Mach numbers on the ascent")
    r.add_argument("case"); r.add_argument("mach", type=float, nargs="+"); r.set_defaults(f=cmd_conditions)

    r = sub.add_parser("static", help="generalised forces of an alpha run minus an alpha-0 run, against piston theory")
    r.add_argument("case"); r.add_argument("run_alpha"); r.add_argument("run_zero"); r.add_argument("-n", type=int, default=4)
    r.set_defaults(f=cmd_static)

    r = sub.add_parser("gaf", help="Q(k) at one Mach from the pulse runs (one per mode)")
    r.add_argument("case"); r.add_argument("runs", nargs="+"); r.add_argument("--zero", help="the A = 0 run")
    r.add_argument("-o", "--out", required=True, help="gaf_M<mach>.npz"); r.add_argument("--plot")
    r.set_defaults(f=cmd_gaf)

    r = sub.add_parser("linearity", help="compare the Q column of the A and 2A pulse runs of one mode")
    r.add_argument("run_a"); r.add_argument("run_b"); r.add_argument("--zero"); r.add_argument("--fmax", type=float, default=1500)
    r.add_argument("--plot"); r.set_defaults(f=cmd_linearity)

    r = sub.add_parser("flutter", help="Roger RFA of the GAF files and the flutter margin over the flight")
    r.add_argument("case"); r.add_argument("gafs", nargs="+"); r.add_argument("--modes", type=int, default=4)
    r.add_argument("--g-file", help="tap.json written by 'tap'"); r.add_argument("-o", "--out", default="results/flutter")
    r.set_defaults(f=cmd_flutter)

    r = sub.add_parser("tap", help="frequencies and loss factors from a tap-test recording (WAV or CSV)")
    r.add_argument("case"); r.add_argument("recording"); r.add_argument("--modes", type=int, default=3)
    r.add_argument("-o", "--out", default="tap.json"); r.add_argument("--plot"); r.set_defaults(f=cmd_tap)

    r = sub.add_parser("lco", help="energy per cycle of the harmonic runs at several amplitudes")
    r.add_argument("case"); r.add_argument("runs", nargs="+"); r.add_argument("--g-file")
    r.add_argument("--q-scale", type=float, default=1.0, help="q_flight / q_cfd of the base run")
    r.add_argument("--plot"); r.set_defaults(f=cmd_lco)

    r = sub.add_parser("sweep", help="forced-run sweep file and its flutterModes files")
    r.add_argument("case"); r.add_argument("base", nargs="+", help="sweep files with the base runs")
    r.add_argument("-o", "--out", required=True); r.add_argument("--modes-dir", default="~/flutter/aconcagua")
    r.add_argument("--modes", type=int, nargs="+", default=[1, 2, 3, 4]); r.add_argument("--amplitude", type=float, default=5e-4)
    r.add_argument("--fmax", type=float, default=1500.0); r.add_argument("--mesh", default="h0")
    r.add_argument("--linearity-machs", type=float, nargs="+", default=[1.05, 1.86])
    r.add_argument("--base-solver", default="pimple", help="solver of the base runs to start from (pimple or central)")
    r.set_defaults(f=cmd_sweep)

    r = sub.add_parser("check", help="check a forced-motion run directory before running it")
    r.add_argument("run"); r.set_defaults(f=cmd_check)

    a = p.parse_args(argv)
    a.argv = argv[1:]
    a.f(a)
