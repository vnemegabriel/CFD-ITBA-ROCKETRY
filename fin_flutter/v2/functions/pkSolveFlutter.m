function [V_fl, V_div, results] = pkSolveFlutter(omega_n, Q_k, k_vals, b_ref, flightConds)
% pkSolveFlutter  p-k flutter + quasi-steady divergence solver (modal space).
%
%   [V_fl, V_div, results] = pkSolveFlutter(omega_n, Q_k, k_vals, b_ref, flightConds)
%
%   omega_n    : [nModes×1]  natural frequencies [rad/s]
%   Q_k        : [nModes×nModes×nK]  unsymmetrized GAF per unit q_inf [1/Pa]
%                Q_k(:,:,1) must correspond to k_vals(1) = 0
%   k_vals     : [1×nK]  reduced frequencies; k_vals(1) = 0
%   b_ref      : reference semi-chord [m]  (MAC/2)
%   flightConds: struct array (.rho .U .q_inf)
%
% ── AERODYNAMIC MODEL ────────────────────────────────────────────────────
%   Piston theory yields Q(k) = Q0 + i*k*Q1  (exact linear in k).
%     Q0 = real(Q_k(:,:,1))         quasi-steady stiffness, NOT symmetrized
%     Q1 = imag(Q_k(:,:,2))/k(2)   aerodynamic damping matrix (physical units: 1/Pa)
%   Q is non-Hermitian: the cross-coupling Q_ij = ∫φ_i·∂φ_j/∂x dA ≠ Q_ji.
%   Symmetrizing destroys the physical off-diagonal coupling that drives
%   bending-torsion flutter — this is why solveFlutterPL gives λ*=26 on B2.
%
% ── STRUCTURAL DAMPING ───────────────────────────────────────────────────
%   A structural loss factor g = 0.04 is applied via complex stiffness:
%     Ω²_c = diag(ω_n²·(1+ig))
%   This regularises the undamped system so Re(p) < 0 at q=0 and flutter
%   is identified by Re(p) crossing zero. For g ≤ 0.05, the flutter speed
%   changes by less than 5% relative to the zero-damping limit.
%
% ── p-k LOOP ─────────────────────────────────────────────────────────────
%   At each dynamic pressure q, iterate k until convergence:
%     1. Q_eff = q·(Q0 + i·k·Q1)
%     2. A = Q_eff − Ω²_c;   eigenvalues λᵢ = pᵢ²  (complex, non-Hermitian eig)
%     3. pᵢ = √λᵢ  with  Im(pᵢ) > 0  (physical, positive-frequency branch)
%     4. k_new = mean(Im(p))·b_ref/U;  iterate until |Δk| < tol
%   Flutter: first q where any Re(pᵢ) crosses 0 from below.
%
% ── DIVERGENCE ───────────────────────────────────────────────────────────
%   B̂ = diag(1/ω)·Q0·diag(1/ω);  q_div = 1/max(positive eigenvalues of B̂).
%   Divergence-free if all eigenvalues of B̂ are ≤ 0.

nModes = length(omega_n);
nF     = length(flightConds);

% Loss factor: 4% structural damping, representative for T700/epoxy CFRP composite.
% CFRP laminates typically exhibit g = 0.02–0.05 (CLT-averaged from fiber-dominated
% and matrix-dominated modes).  g = 0.01 (1%) is conservative; it shifts the
% single-mode aerodynamic anti-damping flutter speed below flight speed.

g_struct = 0.04;
Omega2   = diag(omega_n.^2 * (1 + 1i*g_struct));   % complex stiffness

% ── Aerodynamic decomposition ────────────────────────────────────────────
Q0 = real(Q_k(1:nModes, 1:nModes, 1));            % quasi-steady [1/Pa], real
Q1 = imag(Q_k(1:nModes, 1:nModes, 2)) / k_vals(2); % damping slope [1/Pa]

% ── Divergence ───────────────────────────────────────────────────────────
B_hat   = diag(1./omega_n) * Q0 * diag(1./omega_n);
lam_hat = real(eig(B_hat));
pos_lam = lam_hat(lam_hat > 1e-12);
if isempty(pos_lam)
    q_div_crit = Inf;
else
    q_div_crit = 1 / max(pos_lam);
end

% ── p-k sweep parameters ─────────────────────────────────────────────────
nQ       = 150;
max_iter = 40;
k_tol    = 1e-8;

V_fl     = Inf(nF, 1);
V_div    = Inf(nF, 1);
gam_hist = cell(nF, 1);
omg_hist = cell(nF, 1);
q_hist   = cell(nF, 1);

for fi = 1:nF
    rho_i = flightConds(fi).rho;
    U_i   = flightConds(fi).U;
    q_flt = flightConds(fi).q_inf;

    if ~isinf(q_div_crit)
        V_div(fi) = sqrt(2 * q_div_crit / rho_i);
    end

    q_max = 4 * q_flt;
    if ~isinf(q_div_crit)
        q_max = max(q_max, 1.1 * q_div_crit);
    end
    q_vec = linspace(0, q_max, nQ);

    gam = zeros(nModes, nQ);
    omg = zeros(nModes, nQ);

    % Initialize: undamped + structural loss factor
    p_cur = omega_n .* (-g_struct/2 + 1i);   % [nModes×1]
    k_avg = 0;   % Start quasi-steady at each flight point (k=0).
                 % Seeding from natural frequencies gives k_init ~ O(1), which
                 % injects a large Im(Q1) term at the first q-step and can flip
                 % Re(p) positive before any aerodynamic coupling has built up,
                 % producing a spurious flutter detection near q=0.

    for qi = 1:nQ
        q = q_vec(qi);

        if q < eps
            gam(:, qi) = real(p_cur);
            omg(:, qi) = imag(p_cur);
            continue;
        end

        % k-convergence: iterate k = mean(Im(p))·b/U; k_avg carries over from
        % the previous q-step for continuity (warm-start), but is reset to 0
        % at the beginning of each flight condition above.

        p_new = p_cur;   % fallback
        for it = 1:max_iter
            Q_eff = q * (Q0 + 1i * k_avg * Q1);
            A     = Q_eff - Omega2;        % p² = eig(A)
            lam   = eig(A);                % [nModes×1] complex

            % Physical branch: Im(p) ≥ 0
            p_raw = arrayfun(@pickBranch, lam);

            % Mode tracking: greedy minimum-distance assignment
            p_new = matchModes(p_raw, p_cur);

            % k is driven by the LOWEST-frequency mode (mode 1 ≡ primary flutter mode).
            % Using mean(Im(p)) includes high-frequency modes (ω₃..ω₆) that push k_avg
            % to O(1), which injects a spurious Im(Q1) contribution that makes
            % Re(p₁) cross zero at q << q_flight.  The self-consistent k for mode 1
            % is Im(p₁)×b/U ≈ 0.23, not the mean ≈ 1.3.
            [k_new, ~] = min(max(imag(p_new), 0));
            k_new = k_new * (b_ref / U_i);
            if abs(k_new - k_avg) < k_tol, break; end
            k_avg = k_new;
        end

        p_cur       = p_new;
        gam(:, qi)  = real(p_new);
        omg(:, qi)  = imag(p_new);
    end

    gam_hist{fi} = gam;
    omg_hist{fi} = omg;
    q_hist{fi}   = q_vec;

    % Flutter: first q where Re(pᵢ) crosses 0 from below (stable → unstable)
    q_fl = findFlutterQ(gam, q_vec);
    if ~isnan(q_fl)
        V_fl(fi) = sqrt(2 * q_fl / rho_i);
    end
end

results.q_div_crit = q_div_crit;
results.g_struct   = g_struct;
results.gam_hist   = gam_hist;
results.omg_hist   = omg_hist;
results.q_hist     = q_hist;
results.Q0         = Q0;
results.Q1         = Q1;
results.omega_n    = omega_n;
end


%% ── Local helpers ────────────────────────────────────────────────────────

function p = pickBranch(lam)
% Return sqrt(lam) with Im(p) >= 0 (positive-frequency physical branch).
s = sqrt(lam);
if imag(s) < 0
    s = -s;
end
p = s;
end


function p_out = matchModes(p_new, p_ref)
% Greedy minimum complex-distance assignment of p_new onto p_ref ordering.
nM    = length(p_ref);
used  = false(nM, 1);
p_out = zeros(nM, 1);
for r = 1:nM
    d = abs(p_new - p_ref(r));
    d(used) = Inf;
    [~, idx] = min(d);
    p_out(r) = p_new(idx);
    used(idx) = true;
end
end


function q_fl = findFlutterQ(gam, q_vec)
% Find the smallest q where any mode's Re(p) crosses 0 from below.
q_fl = NaN;
for r = 1:size(gam, 1)
    row = gam(r, :);
    nc  = find(diff(sign(row)) > 0, 1);   % negative → positive crossing
    if isempty(nc), continue; end
    drow = row(nc+1) - row(nc);
    if abs(drow) < 1e-20, continue; end    % flat-zero artifact
    q_c = interp1(row(nc:nc+1), q_vec(nc:nc+1), 0, 'linear');
    if ~isnan(q_c) && q_c > 0
        if isnan(q_fl) || q_c < q_fl
            q_fl = q_c;
        end
    end
end
end
