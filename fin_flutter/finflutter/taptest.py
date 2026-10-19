"""Tap test: natural frequencies and damping of the fin from a recording of free decays."""
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.io import wavfile


def load(path, channel=0):
    """(fs, x) from a WAV file or a CSV with columns time [s], acceleration (any unit)."""
    path = Path(path)
    if path.suffix.lower() == ".wav":
        fs, x = wavfile.read(path)
        x = x.astype(float)
        if x.ndim > 1:
            x = x[:, channel]
        return float(fs), x / np.abs(x).max()
    d = np.loadtxt(path, delimiter=",", comments="#", ndmin=2)
    t, x = d[:, 0], d[:, 1 + channel]
    fs = 1.0 / np.median(np.diff(t))
    return fs, np.interp(np.arange(t[0], t[-1], 1 / fs), t, x)


def taps(x, fs, threshold=0.3, min_gap=0.3):
    """Start indices of taps: the envelope crosses threshold * max, at least min_gap s apart."""
    env = np.abs(signal.hilbert(x))
    env = signal.filtfilt(*signal.butter(2, 200 / (fs / 2)), env)
    above = np.flatnonzero(env > threshold * env.max())
    starts = [above[0]] if len(above) else []
    for i in above[1:]:
        if i - starts[-1] > min_gap * fs:
            starts.append(i)
    return np.array(starts, int)


def peaks(x, fs, fmin=50, fmax=4000, n=4, prominence_db=10.0):
    """The n most prominent spectral peaks of one decay, in Hz. No window: the decay tapers itself."""
    X = 20 * np.log10(np.abs(np.fft.rfft(x, 4 * len(x))) + 1e-12)
    f = np.fft.rfftfreq(4 * len(x), 1 / fs)
    band = (f > fmin) & (f < fmax)
    idx, prop = signal.find_peaks(X[band], prominence=prominence_db)
    order = np.argsort(prop["prominences"])[::-1][:n]
    return np.sort(f[band][idx[order]])


def damping(x, fs, f0, bw=0.12, hi=0.8, lo=0.1):
    """zeta of the mode at f0: band-pass +-bw*f0, Hilbert envelope, slope of ln(env) between hi and lo of its peak."""
    sos = signal.butter(4, [(1 - bw) * f0, (1 + bw) * f0], btype="band", fs=fs, output="sos")
    y = signal.sosfiltfilt(sos, x)
    env = np.abs(signal.hilbert(y))
    i0 = int(np.argmax(env))
    e = env[i0:]
    sel = np.flatnonzero((e < hi * e[0]) & (e > lo * e[0]))
    if len(sel) < 10:
        return np.nan, env[i0]
    sel = sel[: np.argmax(np.diff(sel) > 1) + 1] if np.any(np.diff(sel) > 1) else sel
    t = sel / fs
    slope = np.polyfit(t, np.log(e[sel]), 1)[0]
    return -slope / (2 * np.pi * f0), env[i0]


def analyse(path, n_modes=3, window=0.25, channel=0):
    """Per tap: frequencies and damping of the n_modes strongest peaks. Returns a list of dicts."""
    fs, x = load(path, channel)
    out = []
    for s in taps(x, fs):
        seg = x[s: s + int(window * fs)]
        for f0 in peaks(seg, fs, n=n_modes):
            z, a = damping(seg, fs, f0)
            out.append({"tap_s": s / fs, "f_Hz": f0, "zeta": z, "g": 2 * z, "amplitude": a})
    return out


def summary(rows, f_model, tol=0.15):
    """Group the tap results by the nearest model frequency (within tol); mean and spread per mode."""
    res = []
    for k, fm in enumerate(f_model):
        r = [x for x in rows if abs(x["f_Hz"] - fm) < tol * fm and np.isfinite(x["zeta"])]
        if r:
            f = np.array([x["f_Hz"] for x in r]); g = np.array([x["g"] for x in r]); a = np.array([x["amplitude"] for x in r])
            slope = np.polyfit(a, g, 1)[0] if len(r) > 2 and np.ptp(a) > 0 else np.nan
            res.append({"mode": k + 1, "f_model": fm, "n": len(r), "f_mean": f.mean(), "f_std": f.std(),
                        "g_mean": g.mean(), "g_std": g.std(), "dg_damp": slope})
        else:
            res.append({"mode": k + 1, "f_model": fm, "n": 0})
    return res
