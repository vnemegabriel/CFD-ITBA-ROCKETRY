function [V_fl, V_div, results] = solveFlutterPL(omega_n, Q_k, k_vals, flightConds)
% solveFlutterPL  p-L flutter + static divergence solver (modal-space).
%
%   [V_fl, V_div, results] = solveFlutterPL(omega_n, Q_k, k_vals, flightConds)
%
%   omega_n      : [nModes × 1]  natural frequencies [rad/s]
%   Q_k          : [nModes × nModes × nK]  GAF normalised per unit q [1/Pa]
%                  Q_k(:,:,1)  MUST correspond to k_vals(1) = 0  (quasi-steady)
%   k_vals       : [1 × nK]  reduced frequencies; k_vals(1) = 0 required
%   flightConds  : struct array  (.Mach .a .rho .U .q_inf)
%
%   V_fl   : [nF×1]  flutter    speed [m/s]  (Inf = stable within sweep)
%   V_div  : [nF×1]  divergence speed [m/s]  (Inf = divergence-free design)
%   results: struct  (flutter_margin, div_margin, eig_history_flutter, ...)
%
% ── DIVERGENCE (static, k = 0) ──────────────────────────────────────────
%   Uses Q0 = real(Q_k(:,:,1)) — quasi-steady (k=0) aerodynamic stiffness.
%
%   Divergence condition:   det( K_modal − q_div · Q0 ) = 0
%
%   Transformed to a standard symmetric eigenvalue to avoid MATLAB's
%   generalised-eig instability when Q0 is indefinite:
%
%       B̂ = diag(1/ω) · Q0 · diag(1/ω)          [nModes × nModes]
%       eig(B̂) = λ̂ᵢ  →  q_div,i = 1/λ̂ᵢ   for λ̂ᵢ > 0
%
%   Critical divergence: q_div = 1 / max( λ̂ᵢ  [positive] )
%   If all λ̂ᵢ ≤ 0 (wash-out design): V_div = Inf (divergence-free).
%
% ── FLUTTER (dynamic, k > 0) ────────────────────────────────────────────
%   Uses Q_dyn = real(mean(Q_k(:,:,2:end), 3)) — average over k > 0.
%   Sweeps q from 0 to 4×q_flight; detects first sign change of min eig of
%   K_ae(q) = K_modal − q · Q_dyn.
%
% ── STABILITY MARGINS ───────────────────────────────────────────────────
%   flutter_margin = q_flutter_crit / q_flight   [> 1 = safe]
%   div_margin     = q_div_crit     / q_flight   [> 1 = safe, Inf = free]

% -------------------------------------------------------------------------
nModes = length(omega_n);
nF     = length(flightConds);
nK     = size(Q_k, 3);

K_modal = diag(omega_n.^2);          % [nModes × nModes]

% =========================================================================
% DIVERGENCE — quasi-steady k = 0 analysis
% =========================================================================
Q0     = real(Q_k(1:nModes, 1:nModes, 1));   % k = 0 slice (purely real)
Q0_sym = (Q0 + Q0') / 2;

% Standard symmetric eigenvalue for divergence:
%   Derived from K_modal*φ = q_div*Q0*φ  →  B̂*ψ = (1/q_div)*ψ
%   where B̂ = diag(1/ω)*Q0*diag(1/ω)
B_hat   = diag(1./omega_n) * Q0_sym * diag(1./omega_n);
B_hat   = (B_hat + B_hat') / 2;
lam_hat = real(eig(B_hat));           % [nModes × 1]

pos_lam = lam_hat(lam_hat > 1e-12);   % positive → finite divergence q

if isempty(pos_lam)
    q_div_crit = Inf;    % all Q0 eigenvalues ≤ 0: divergence-free (wash-out)
else
    q_div_crit = 1 / max(pos_lam);    % smallest critical q = 1/largest λ̂
end

% =========================================================================
% FLUTTER — dynamic k > 0 analysis
% =========================================================================
if nK > 1
    % Retain the full complex average — imaginary part carries aerodynamic damping.
    % K_ae = K_modal - q*Q_dyn remains Hermitian (K_modal real-symmetric,
    % Q_dyn complex Hermitian) so eigenvalues are real.  Discarding Im(Q_dyn)
    % loses the frequency-dependent stiffness shift from unsteady loading.
    % NOTE: this is still a quasi-steady check (no k–ω iteration); classical
    % bending-torsion flutter via mode coalescence requires a p-k loop —
    % see loewnerInterpolation.m for the intended extension.
    Q_dyn = mean(Q_k(1:nModes, 1:nModes, 2:end), 3);
else
    Q_dyn = Q0;
end
Q_dyn = (Q_dyn + Q_dyn') / 2;  % enforce Hermitian (guard against floating-point asymmetry)

% Analytical flutter margin (extrapolation when no crossing found in sweep)
lam_Qdyn     = real(eig(Q_dyn));
lam_Qdyn_max = max(lam_Qdyn);
if lam_Qdyn_max > 1e-12
    q_flutter_crit = min(omega_n.^2) / lam_Qdyn_max;
else
    q_flutter_crit = Inf;   % flutter-free (wash-out)
end

% =========================================================================
% Allocate outputs
% =========================================================================
V_fl           = Inf(nF, 1);
V_div          = Inf(nF, 1);
flutter_margin = NaN(nF, 1);
div_margin     = NaN(nF, 1);

results.eig_history_flutter = cell(nF, 1);
results.q_sweep             = cell(nF, 1);
results.omega_n             = omega_n;
results.Q0_sym              = Q0_sym;
results.Q_dyn               = Q_dyn;
results.q_div_crit          = q_div_crit;
results.q_flutter_crit      = q_flutter_crit;
results.lam_hat_diverge     = lam_hat;    % eigenvalues of B̂  (= 1/q_div per mode)
results.lam_Qdyn            = lam_Qdyn;  % eigenvalues of Q_dyn

nQ = 80;

% =========================================================================
% Per-flight-point loop
% =========================================================================
for fi = 1:nF
    rho_i = flightConds(fi).rho;
    q_fl  = flightConds(fi).q_inf;

    % ── Divergence ───────────────────────────────────────────────────────
    if isinf(q_div_crit)
        V_div(fi)      = Inf;
        div_margin(fi) = Inf;
    else
        V_div(fi)      = sqrt(2 * q_div_crit / rho_i);
        div_margin(fi) = q_div_crit / q_fl;
    end

    % ── Flutter: q-sweep on Q_dyn ────────────────────────────────────────
    q_max_sweep = 4 * q_fl;
    if ~isinf(q_div_crit)
        q_max_sweep = max(q_max_sweep, q_div_crit * 1.05);
    end
    q_vec = linspace(0, q_max_sweep, nQ);

    eig_lam = zeros(nModes, nQ);
    for qi = 1:nQ
        K_ae = K_modal - q_vec(qi) * Q_dyn;
        eig_lam(:, qi) = sort(real(eig(K_ae)));
    end
    results.eig_history_flutter{fi} = eig_lam;
    results.q_sweep{fi}             = q_vec;

    fl_found = false;
    for m = 1:nModes
        row = eig_lam(m, :);
        nc  = find(diff(sign(row)) < 0, 1, 'first');
        if isempty(nc), continue; end
        q_found = interp1(row(nc:nc+1), q_vec(nc:nc+1), 0, 'linear');
        if ~isnan(q_found) && q_found > 0
            V_fl(fi)           = sqrt(2 * q_found / rho_i);
            flutter_margin(fi) = q_found / q_fl;
            fl_found = true;
            break;
        end
    end
    if ~fl_found
        flutter_margin(fi) = q_flutter_crit / q_fl;
    end
end

% =========================================================================
% Pack results
% =========================================================================
results.V_fl             = V_fl;
results.V_div            = V_div;
results.flutter_margin   = flutter_margin;
results.div_margin       = div_margin;
results.stability_margin = flutter_margin;   % legacy alias
end
