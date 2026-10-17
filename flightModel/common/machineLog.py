#!/usr/bin/env python3
"""State of the machine during a run, so cost numbers can be compared.

    machineLog.py snapshot                       static state, as foamDictionary entries
    machineLog.py sample <csv> [seconds]         one CSV line every interval until SIGTERM
    machineLog.py summary <csv> [perf] [energy]  totals, as foamDictionary entries

perf and energy are the -x, -I outputs of perf stat.  See docs/MEASUREMENTS.md.
"""
import csv
import glob
import os
import platform
import signal
import subprocess
import sys
import time


def read(path, default=""):
    try:
        with open(path) as f:
            return f.read().strip()
    except OSError:
        return default


def hwmon():
    """name -> list of hwmon directories (spd5118 appears once per DIMM)."""
    out = {}
    for d in sorted(glob.glob("/sys/class/hwmon/hwmon*")):
        out.setdefault(read(d + "/name"), []).append(d)
    return out


def milli(path):
    v = read(path)
    return "" if v == "" else "%.1f" % (int(v) / 1000)


def cpu_mhz():
    vals = []
    for c in glob.glob("/sys/devices/system/cpu/cpu[0-9]*/cpufreq"):
        v = read(c + "/cpuinfo_avg_freq") or read(c + "/scaling_cur_freq")
        if v:
            vals.append(int(v) / 1000)
    return vals


def fan_level():
    for line in read("/proc/acpi/ibm/fan").splitlines():
        if line.startswith("level:"):
            return line.split()[1]
    return "unknown"


def snapshot():
    cpu0 = "/sys/devices/system/cpu/cpu0/cpufreq/"
    lscpu = subprocess.run(["lscpu"], capture_output=True, text=True).stdout
    info = dict(l.split(":", 1) for l in lscpu.splitlines() if ":" in l)
    mem = [l for l in read("/proc/meminfo").splitlines() if l.startswith("MemTotal")]
    try:
        mpi = subprocess.run(["ompi_info", "--version"], capture_output=True,
                             text=True).stdout.splitlines()[0]
    except (OSError, IndexError):
        mpi = "unknown"
    wm = " ".join(os.environ.get(k, "") for k in (
        "WM_PROJECT_VERSION", "WM_COMPILER", "WM_COMPILE_OPTION",
        "WM_PRECISION_OPTION", "WM_LABEL_SIZE", "WM_MPLIB"))
    e = {
        "kernel": platform.release(),
        "sockets": info.get("Socket(s)", "").strip(),
        "coresPerSocket": info.get("Core(s) per socket", "").strip(),
        "threadsPerCore": info.get("Thread(s) per core", "").strip(),
        "memTotalkB": mem[0].split()[1] if mem else "",
        "governor": read(cpu0 + "scaling_governor"),
        "epp": read(cpu0 + "energy_performance_preference"),
        "boost": read(cpu0 + "boost", read("/sys/devices/system/cpu/cpufreq/boost")),
        "scalingMaxMHz": int(read(cpu0 + "scaling_max_freq", "0")) // 1000,
        "platformProfile": read("/sys/firmware/acpi/platform_profile", "none"),
        "fanLevel": fan_level(),
        "perfParanoid": read("/proc/sys/kernel/perf_event_paranoid"),
        "openfoam": wm,
        "mpi": mpi,
    }
    for k, v in e.items():
        print('    %-15s "%s";' % (k, v))


def sample(out, interval):
    mon = hwmon()
    dimms = [d + "/temp1_input" for d in mon.get("spd5118", [])]
    tctl = [d + "/temp1_input" for d in mon.get("k10temp", []) + mon.get("coretemp", [])][:1]
    nvme = [d + "/temp1_input" for d in mon.get("nvme", [])][:1]
    fans = sorted(glob.glob(mon.get("thinkpad", [""])[0] + "/fan*_input"))
    head = (["wall", "monotonic", "boottime", "mhzMean", "mhzMin", "mhzMax", "load1",
             "tCpu"] + ["tDimm%d" % i for i in range(len(dimms))] + ["tNvme"]
            + ["fan%d" % i for i in range(len(fans))] + ["profile"])
    stop = []
    signal.signal(signal.SIGTERM, lambda *a: stop.append(1))
    signal.signal(signal.SIGINT, lambda *a: stop.append(1))
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(head)
        while True:
            mhz = cpu_mhz() or [0]
            w.writerow(["%.0f" % time.time(),
                        "%.1f" % time.clock_gettime(time.CLOCK_MONOTONIC),
                        "%.1f" % time.clock_gettime(time.CLOCK_BOOTTIME),
                        "%.0f" % (sum(mhz) / len(mhz)), "%.0f" % min(mhz), "%.0f" % max(mhz),
                        "%.2f" % os.getloadavg()[0]]
                       + [milli(p) for p in tctl + dimms + nvme]
                       + [read(p) for p in fans]
                       + [read("/sys/firmware/acpi/platform_profile", "none")])
            f.flush()
            if stop:
                return
            t = time.monotonic() + interval
            while not stop and time.monotonic() < t:
                time.sleep(1)


def perf_totals(path):
    """Sum the per-interval counts of a perf stat -x, -I file."""
    tot = {}
    for row in csv.reader(l for l in open(path) if l.strip() and not l.startswith("#")):
        if len(row) < 4:
            continue
        try:
            tot[row[3]] = tot.get(row[3], 0.0) + float(row[1])
        except ValueError:
            pass
    return tot


def summary(machine, perf=None, energy=None):
    e = {}
    rows = list(csv.DictReader(open(machine)))
    if rows:
        a, b = rows[0], rows[-1]
        e["samples"] = len(rows)
        e["suspendedSeconds"] = "%.0f" % max(0.0, (
            (float(b["boottime"]) - float(a["boottime"]))
            - (float(b["monotonic"]) - float(a["monotonic"]))))
        mhz = [float(r["mhzMean"]) for r in rows]
        e["mhzMean"] = "%.0f" % (sum(mhz) / len(mhz))
        e["mhzMinOfMean"] = "%.0f" % min(mhz)
        for k in [c for c in rows[0] if c.startswith("t")]:
            v = [float(r[k]) for r in rows if r[k] != ""]
            if v:
                e[k + "Max"] = "%.1f" % max(v)
        e["profiles"] = '"%s"' % " ".join(sorted({r["profile"] for r in rows}))
    if perf and os.path.exists(perf):
        t = perf_totals(perf)
        tc = t.get("task-clock", 0) * 1e-3
        cyc = t.get("cycles", 0)
        ins = t.get("instructions", 0)
        if tc > 0 and cyc > 0:
            e["cpuSeconds"] = "%.1f" % tc
            e["cycles"] = "%.6g" % cyc
            e["instructions"] = "%.6g" % ins
            e["GHzEffective"] = "%.3f" % (cyc / tc / 1e9)
            e["IPC"] = "%.3f" % (ins / cyc)
            if t.get("stalled-cycles-frontend"):
                e["frontendStallFraction"] = "%.3f" % (t["stalled-cycles-frontend"] / cyc)
            if t.get("cache-misses"):
                e["cacheMisses"] = "%.6g" % t["cache-misses"]
            e["cpuMigrations"] = "%.0f" % t.get("cpu-migrations", 0)
            e["contextSwitches"] = "%.0f" % t.get("context-switches", 0)
    if energy and os.path.exists(energy):
        j = sum(v for k, v in perf_totals(energy).items() if "energy-pkg" in k)
        if j > 0:
            e["energyPkgJ"] = "%.0f" % j
    for k, v in e.items():
        print("    %-21s %s;" % (k, v))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "snapshot":
        snapshot()
    elif cmd == "sample":
        sample(args[0], float(args[1]) if len(args) > 1 else 30)
    elif cmd == "summary":
        summary(*args)
    else:
        sys.exit(__doc__)
