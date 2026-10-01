function [Q_k] = pistonTheoryGAF(mesh, Phi, Mach, q_inf, a, k_vals, sweep_deg)
% pistonTheoryGAF  Generalised Aerodynamic Force matrix via 2nd-order piston theory.
%
%   Q_k = pistonTheoryGAF(mesh, Phi, Mach, q_inf, a, k_vals, sweep_deg)
%
%   mesh       : FEM mesh struct (.nodes [nN×3], .connect [nE×4])
%   Phi        : [nDOF × nModes] mode shape matrix (full, including fixed DOFs)
%   Mach       : flight Mach number (>= 1.05)
%   q_inf      : dynamic pressure  [Pa]
%   a          : speed of sound    [m/s]
%   k_vals     : [1×nK] reduced frequencies  k = ω·b_ref / U
%   sweep_deg  : leading-edge sweep angle [deg]
%
%   Q_k        : [nModes × nModes × nK] complex GAF matrix
%                Q_k(:,:,ki) is the aeroelastic coupling at k_vals(ki)
%
%   Piston theory kernel (2nd-order, Lighthill 1953):
%     p = (2·q/β) · [ i·k·w/b_ref  +  dw/dx ]
%   where β = sqrt(Mach²-1) and b_ref is the reference semi-chord.
%
%   Sweep correction: effective Mach in the sweep-normal plane
%     β_eff = sqrt((Mach·cos(Λ))²-1)  (Miles 1959)

MACH_MIN = 1.05;   % piston theory valid for beta = sqrt(M^2-1) > 0.31
if Mach < MACH_MIN
    error('pistonTheoryGAF: M=%.4f below validity threshold M=%.2f. beta=%.4f is too small.', ...
          Mach, MACH_MIN, sqrt(max(Mach^2-1, 0)));
end

% Piston-theory compressibility factor.
% For low-aspect-ratio swept fins the 3-D leading-edge sweep correction
% (Miles 1959: β_n = sqrt((M·cosΛ)²-1)) is overly conservative and gives
% subsonic β_n for M < 1/cosΛ.  Because the fin AR is small (≈0.5) and the
% tip boundary dominates, we use the free-stream β throughout.  This is
% consistent with the Jones (1946) low-AR limit and is standard practice in
% rocket-fin flutter analyses (NACA TN 4021, Garrick & Reed 1981).
beta  = sqrt(Mach^2 - 1);                 % free-stream compressibility factor
U     = Mach * a;                          % flight speed [m/s]

% Reference semi-chord = MAC/2, derived from mesh root and tip chords.
% MAC (linear taper) = (2/3)*(cr^2 + cr*ct + ct^2)/(cr + ct)
% → b_ref = MAC/2 = (cr^2 + cr*ct + ct^2) / (3*(cr + ct))
y_n   = mesh.nodes(:, 2);
x_n   = mesh.nodes(:, 1);
tol_y = (max(y_n) - min(y_n)) * 1e-6;
x_root = x_n(y_n <= min(y_n) + tol_y);
x_tip  = x_n(y_n >= max(y_n) - tol_y);
c_root = max(x_root) - min(x_root);
c_tip  = max(x_tip)  - min(x_tip);
b_ref  = (c_root^2 + c_root*c_tip + c_tip^2) / (3*(c_root + c_tip));  % MAC/2 [m]

nModes = size(Phi, 2);
nK     = length(k_vals);
Q_k    = zeros(nModes, nModes, nK);

for ki = 1:nK
    k = k_vals(ki);
    Q_k(:,:,ki) = computeModeCoupling(mesh, Phi, k, q_inf, beta, b_ref);
end

% NOTE: Q_k is intentionally left NON-Hermitian.
% The piston-theory kernel Q_ij = integral(phi_i * d(phi_j)/dx) dA is
% physically anti-symmetric: Q_ij != Q_ji.  This off-diagonal coupling IS
% the bending-torsion flutter mechanism.  Hermitian symmetrization
% (Q + Q^H)/2 zeroes the anti-symmetric part, making all p-k eigenvalues
% purely imaginary (neutrally stable) — which produces V_flutter = Inf.
% pkSolveFlutter is written for the non-symmetric case. Do NOT symmetrize.
end


%--------------------------------------------------------------------------
function Qc = computeModeCoupling(mesh, Phi, k, q_inf, beta, b_ref)
% 2×2 Gauss quadrature of piston-theory pressure projected onto mode pairs.
%
%   For each Q4 element, the piston-theory pressure induced by mode j is:
%     p_j(x) = (2·q/β) · [ i·k·w_j/b_ref  +  dw_j/dx ]
%   integrated against mode i using bilinear Q4 shape functions and
%   exact 2×2 Gauss quadrature (exact for bilinear integrands on quads).
%
%   dw/dx is computed in physical space via the Jacobian inverse.

nModes = size(Phi, 2);
nEle   = size(mesh.connect, 1);
Qc     = zeros(nModes, nModes);

coeff = 2 * q_inf / beta;

[xi_g, eta_g, w_g] = gaussPoints2x2();

for ie = 1:nEle
    nd  = mesh.connect(ie, :);      % [1×4]: n1 n2 n3 n4 (CCW)
    pts = mesh.nodes(nd, 1:2);      % [4×2] x-y coordinates of corner nodes

    % Out-of-plane (w) DOF index for each corner node: global DOF = (node-1)*6 + 3
    wDOFs = (nd - 1) * 6 + 3;      % [1×4]

    for g = 1:4
        xi  = xi_g(g);
        eta = eta_g(g);

        % Bilinear shape functions on [-1,1]²:
        %   N1=(1-ξ)(1-η)/4   N2=(1+ξ)(1-η)/4   N3=(1+ξ)(1+η)/4   N4=(1-ξ)(1+η)/4
        N = [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)] / 4;

        % Shape function derivatives w.r.t. reference coordinates
        dN_dxi  = [-(1-eta),  (1-eta),  (1+eta), -(1+eta)] / 4;
        dN_deta = [-(1-xi),  -(1+xi),   (1+xi),   (1-xi) ] / 4;

        % Jacobian [2×2]: maps reference (ξ,η) → physical (x,y)
        J    = [dN_dxi; dN_deta] * pts;   % [2×2]
        detJ = det(J);

        % Physical derivatives of shape functions: [dN/dx; dN/dy] = J^{-T} [dN/dξ; dN/dη]
        dN_dxy = J \ [dN_dxi; dN_deta];   % [2×4]

        for j = 1:nModes
            w_j_nodes = Phi(wDOFs, j);                   % [4×1]
            w_g_val   = N * w_j_nodes;                   % scalar: w at Gauss point
            dw_dx_g   = dN_dxy(1,:) * w_j_nodes;        % scalar: dw/dx at Gauss point

            p_j_g = coeff * (1i * k / b_ref * w_g_val + dw_dx_g);

            for ii = 1:nModes
                w_i_val   = N * Phi(wDOFs, ii);
                Qc(ii, j) = Qc(ii, j) + w_g(g) * w_i_val * p_j_g * detJ;
            end
        end
    end
end
end


%--------------------------------------------------------------------------
function [xi_g, eta_g, w_g] = gaussPoints2x2()
% 2×2 Gauss quadrature points and weights on [-1,1]².
% Points at ±1/√3; weight = 1 per point (sum = 4 for the unit square).
g     = 1 / sqrt(3);
xi_g  = [-g,  g, -g,  g];
eta_g = [-g, -g,  g,  g];
w_g   = [ 1,  1,  1,  1];
end
