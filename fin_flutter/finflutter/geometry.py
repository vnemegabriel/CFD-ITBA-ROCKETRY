"""Fin geometry read from flightModel/mesh/makeRocket.py, so the FE fin is the CFD fin."""
import ast
import math
import re
from dataclasses import dataclass
from pathlib import Path

MAKEROCKET = Path(__file__).parents[2] / "flightModel/mesh/makeRocket.py"


def rocket_dims(makerocket=MAKEROCKET, dims=None, **override):
    """Defaults of makeRocket.py's argparse options, then a dims file, then keyword overrides."""
    vals = {}
    for node in ast.walk(ast.parse(Path(makerocket).read_text())):
        if isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "add_argument" and node.args:
            kw = {k.arg: k.value for k in node.keywords}
            if "default" in kw and isinstance(node.args[0], ast.Constant) and node.args[0].value.startswith("--"):
                vals[node.args[0].value[2:]] = eval(compile(ast.Expression(kw["default"]), "", "eval"), {"__builtins__": {}})
    if dims:
        for line in Path(dims).read_text().splitlines():
            line = re.split(r"#|//", line)[0].strip().rstrip(";").strip()
            if line:
                name, value = line.split()
                vals[name] = type(vals[name])(value)
    vals.update(override)
    return vals


@dataclass
class Fin:
    cr: float
    ct: float
    span: float
    sweep: float        # tip LE aft of root LE [m]
    t: float            # max thickness, mid-chord
    rho: float          # LE/TE/tip round radius
    R: float            # body radius = root of the fin
    x_root_le: float    # absolute x of the root LE, nose tip at x = 0
    roll_deg: float     # angle of fin 0 about +x, from +y towards +z

    @property
    def sweep_deg(self):
        return math.degrees(math.atan2(self.sweep, self.span))

    def chord(self, eta):
        return self.cr + eta * (self.ct - self.cr)

    def x_le(self, eta):
        return eta * self.sweep

    def half_thickness(self, x, c):
        """Upper surface z(x) of makeRocket.section(c): LE round, tangent segment, circular arc, mirrored aft."""
        t, rho = self.t, self.rho
        x = min(max(x, 0.0), c)
        if x > c / 2:
            x = c - x
        Rc = (c * c + t * t) / (4 * t)
        d = Rc - t / 2
        Dx, Dz = rho - c / 2, d
        beta = math.atan2(Dz, Dx) - math.acos((Rc - rho) / math.hypot(Dx, Dz))
        T1 = (rho + rho * math.cos(beta), rho * math.sin(beta))
        T2 = (c / 2 + Rc * math.cos(beta), -d + Rc * math.sin(beta))
        if x <= T1[0]:
            return math.sqrt(max(rho * rho - (x - rho) ** 2, 0.0))
        if x <= T2[0]:
            return T1[1] + (x - T1[0]) * (T2[1] - T1[1]) / (T2[0] - T1[0])
        return -d + math.sqrt(Rc * Rc - (x - c / 2) ** 2)

    def thickness(self, x, y):
        """Local thickness at fin-plane point (x from root LE, y from root along the span)."""
        eta = min(max(y / self.span, 0.0), 1.0)
        return 2 * self.half_thickness(x - self.x_le(eta), self.chord(eta))


def fin_from_flightmodel(makerocket=MAKEROCKET, dims=None, **override):
    a = rocket_dims(makerocket, dims, **override)
    x_base = a["noseLD"] * a["diameter"] + a["bodyLength"] + a["tailLength"]
    return Fin(cr=a["finRootChord"], ct=a["finTipChord"], span=a["finSpan"], sweep=a["finSweep"],
               t=a["finThickness"], rho=a["finEdgeRadius"], R=a["diameter"] / 2,
               x_root_le=x_base - a["finFromBase"] - a["finRootChord"], roll_deg=a["finRoll"])
