function [Phi, omega_n] = modalAnalysis(K, M, nModes)
% modalAnalysis  Solve generalized eigenvalue problem (K - ω²·M)·Φ = 0.
%
%   [Phi, omega_n] = modalAnalysis(K, M, nModes)
%
%   K, M    : symmetric positive-semi-definite sparse matrices (reduced system)
%   nModes  : number of modes to extract
%
%   Phi     : [nDOF × nModes] mass-normalised mode shapes
%   omega_n : [nModes × 1] natural frequencies [rad/s], sorted ascending
%
%   Method:
%   - Small systems (nDOF < 600): direct LAPACK full eig — robust against K noise.
%   - Large systems: eigs with a small positive shift (K_reg = K + σ·I) so that
%     the single spurious near-zero eigenvalue from floating-point K drift gets
%     pushed to λ ≈ σ, well below the first physical mode. Request nModes+2 modes
%     and discard any below the physical gap.

nDOF = size(K, 1);

if nDOF <= 2500
    %----------------------------------------------------------------------
    % Direct solver: full generalized eig via LAPACK
    %----------------------------------------------------------------------
    [Phi_all, Omega2_all] = eig(full(K), full(M));
    lam_all = real(diag(Omega2_all));

    % Filter: discard spurious eigenvalues. The spurious mode from K floating-
    % point drift satisfies |λ_spurious| ≪ 1e-10 × max(λ_physical).
    lam_max  = max(lam_all);
    tol_lam  = 1e-10 * lam_max;
    phys_mask = lam_all > tol_lam;

    lam_phys = lam_all(phys_mask);
    Phi_phys = Phi_all(:, phys_mask);
else
    %----------------------------------------------------------------------
    % Iterative solver: eigs with small positive shift
    %
    % K_red may have one spurious near-zero eigenvalue from machine precision.
    % Shift: σ = 1 rad²/s² (below any structural mode but above spurious ≈0).
    % Request nModes+3 to have margin after filtering.
    %----------------------------------------------------------------------
    sigma     = 1.0;                               % [rad²/s²] — corresponds to 0.16 Hz
    % Correct shift for generalised EVP: K·φ = λ·M·φ
    % Adding σ·M (not σ·I) ensures λ_true = λ_shifted - σ exactly.
    % Using σ·I instead would apply a different effective shift to translational
    % vs rotational DOFs (since M is not a multiple of I), which inflates
    % frequencies when rotational mass entries m_rot ≪ m_trans.
    K_reg     = K + sigma * M;

    n_request = min(nModes + 3, nDOF - 1);

    opts_e.disp   = 0;
    opts_e.tol    = 1e-10;
    opts_e.maxit  = 600;

    try
        [Phi_all, Omega2_all] = eigs(K_reg, M, n_request, 'sm', opts_e);
        lam_all = real(diag(Omega2_all)) - sigma;  % undo shift: λ_true = λ_shifted - σ

        % Filter: keep eigenvalues that correspond to physical modes (λ > σ/2)
        % The spurious mode lands at ≈ σ after un-shifting; physical modes are >> σ.
        phys_mask = lam_all > sigma * 0.5;
        lam_phys  = lam_all(phys_mask);
        Phi_phys  = Phi_all(:, phys_mask);

    catch ME
        warning('modalAnalysis:eigsFailed', ...
            'eigs failed (%s). Falling back to full eig.', ME.message);
        [Phi_all, Omega2_all] = eig(full(K_reg), full(M));
        lam_all  = real(diag(Omega2_all)) - sigma;
        phys_mask = lam_all > sigma * 0.5;
        lam_phys  = lam_all(phys_mask);
        Phi_phys  = Phi_all(:, phys_mask);
    end
end

if numel(lam_phys) < nModes
    warning('modalAnalysis:fewModes', ...
        'Only %d physical modes found (requested %d). Check BCs.', ...
        numel(lam_phys), nModes);
    nModes = numel(lam_phys);
end

% Sort ascending
[lam_phys, idx] = sort(lam_phys);
Phi_phys = Phi_phys(:, idx);

% Extract requested modes
omega_n = sqrt(max(lam_phys(1:nModes), 0));
Phi     = Phi_phys(:, 1:nModes);

% Mass-normalise: Φ'·M·Φ = I
for i = 1:nModes
    scale = sqrt(Phi(:,i)' * M * Phi(:,i));
    if scale > 0
        Phi(:,i) = Phi(:,i) / scale;
    end
end
end
