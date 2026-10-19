/*---------------------------------------------------------------------------*\
    Structural modes of fin 0 and the imposed modal signal, read from
    constant/flutterModes (written by fin_flutter: python -m finflutter modes).
    Used by 0/pointDisplacement (case-forced) and FOs/FOgenForces.
\*---------------------------------------------------------------------------*/
#ifndef flutterModes_H
#define flutterModes_H

#include "IFstream.H"
#include "dictionary.H"
#include "Time.H"
#include "unitConversion.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace flutter
{

class Modes
{
public:
    scalar xRootLE, rBody, span, cr, ct, sweep, halfWidth;
    label nXi, nEta;
    vector er, en;
    List<scalarField> w;
    word signal;
    scalarField coef, phase;
    scalar t0, tau, omega, tRamp;

    explicit Modes(const fileName& f)
    {
        IFstream is(f);
        if (!is.good())
        {
            FatalErrorInFunction << "cannot read " << f << exit(FatalError);
        }
        const dictionary d(is);
        xRootLE = d.get<scalar>("xRootLE");
        rBody = d.get<scalar>("rBody");
        span = d.get<scalar>("span");
        cr = d.get<scalar>("rootChord");
        ct = d.get<scalar>("tipChord");
        sweep = d.get<scalar>("sweep");
        halfWidth = d.get<scalar>("halfWidth");
        nXi = d.get<label>("nXi");
        nEta = d.get<label>("nEta");
        const scalar roll = degToRad(d.get<scalar>("rollDeg"));
        er = vector(0, Foam::cos(roll), Foam::sin(roll));
        en = vector(0, -Foam::sin(roll), Foam::cos(roll));
        w = d.get<List<scalarField>>("modes");
        const dictionary& s = d.subDict("signal");
        signal = s.get<word>("type");
        coef = s.get<scalarField>("coef");
        phase = s.getOrDefault<scalarField>("phase", scalarField(coef.size(), Zero));
        t0 = s.getOrDefault<scalar>("t0", 0);
        tau = s.getOrDefault<scalar>("tau", 1);
        omega = s.getOrDefault<scalar>("omega", 0);
        tRamp = s.getOrDefault<scalar>("tRamp", 0);
        if (coef.size() != w.size() || phase.size() != w.size())
        {
            FatalErrorInFunction << f << ": signal coef/phase need one value per mode"
                << exit(FatalError);
        }
    }

    // w of mode k [m per unit modal coordinate] at point p; 0 away from fin 0
    scalar shape(const label k, const vector& p) const
    {
        const vector r(0, p.y(), p.z());
        const scalar s = (r & er) - rBody;
        if (s < -halfWidth || mag(r & en) > halfWidth)
        {
            return 0;
        }
        const scalar eta = min(max(s/span, scalar(0)), scalar(1));
        const scalar c = cr + eta*(ct - cr);
        const scalar xi = min(max((p.x() - xRootLE - eta*sweep)/c, scalar(0)), scalar(1));
        const scalar fi = xi*(nXi - 1), fj = eta*(nEta - 1);
        const label i = min(label(fi), nXi - 2), j = min(label(fj), nEta - 2);
        const scalar a = fi - i, b = fj - j;
        const scalarField& W = w[k];
        return (1 - a)*(1 - b)*W[j*nXi + i] + a*(1 - b)*W[j*nXi + i + 1]
             + (1 - a)*b*W[(j + 1)*nXi + i] + a*b*W[(j + 1)*nXi + i + 1];
    }

    // imposed modal coordinate of mode k at time t
    scalar eta(const label k, const scalar t) const
    {
        if (signal == "pulse")
        {
            return coef[k]*Foam::exp(-sqr((t - t0)/tau));
        }
        if (signal == "harmonic")
        {
            const scalar ramp = t < tRamp
                ? 0.5*(1 - Foam::cos(constant::mathematical::pi*t/tRamp)) : 1;
            return coef[k]*ramp*Foam::sin(omega*t + phase[k]);
        }
        return 0;
    }
};


inline const Modes& modes(const Time& runTime)
{
    static autoPtr<Modes> m;
    if (!m)
    {
        m.reset(new Modes(runTime.globalPath()/"constant"/"flutterModes"));
    }
    return *m;
}

} // End namespace flutter
} // End namespace Foam

#endif
