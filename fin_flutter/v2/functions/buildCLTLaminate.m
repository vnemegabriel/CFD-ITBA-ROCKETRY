function [D_flex, t_total, rho_lam, D_info] = buildCLTLaminate(beta_deg, Vf, matprops)
% buildCLTLaminate  AR1 T700/Epoxy bending stiffness via CLT + Halpin-Tsai.
%
%   [D_flex, t_total, rho_lam, D_info] = buildCLTLaminate(beta_deg, Vf)
%   [D_flex, t_total, rho_lam, D_info] = buildCLTLaminate(beta_deg, Vf, matprops)
%
%   Inputs
%     beta_deg  : global ply rotation for aeroelastic tailoring [deg]   (default 20)
%     Vf        : fibre volume fraction                                   (default 0.50)
%     matprops  : (optional) struct to override nominal UD properties
%                   fields: .E1  .E2  .G12  .nu12  [Pa, Pa, Pa, -]
%
%   Outputs
%     D_flex  : [3×3] CLT bending stiffness matrix [N·m]
%               layout [D11 D12 D16; D12 D22 D26; D16 D26 D66]
%     t_total : total laminate thickness [m]
%     rho_lam : laminate density [kg/m³]  (rule of mixtures)
%     D_info  : struct — individual D terms, ply geometry, material data
%
%   AR1 layup (symmetric, 9 physical groups per half-stack, outer → mid-plane):
%     [DB300±45, GA90R, GA90R, DB300±45, GA90R, DB300±45, GA90R, GA90R, DB300±45]
%
%     DB300 biaxial NCF  → 2 UD sublayers: [+45+β, −45+β]  each of thickness t_DB
%     GA90R woven 0/90   → 1 sublayer:     [β]               of thickness t_GA
%     Full laminate = top half + mirror of top half → symmetric (B-matrix = 0).
%
%   Fabric areal weights (FAW):
%     CARBONODB300  ±45 NCF biaxial  FAW = 300 g/m²
%     CARBONOGA90R  woven 0/90       FAW = 302 g/m²
%     t_DB = FAW_DB/2 / (rho_f · Vf)   [half-FAW per UD sublayer]
%     t_GA = FAW_GA   / (rho_f · Vf)
%
%   Micromechanics (Jones 1999 §3):
%     E1   = Ef1·Vf + Em·(1−Vf)                  [ROM, fibre-dominated]
%     E2   = Em·(1+2η·Vf)/(1−η·Vf)               [Halpin-Tsai, ξ=2]
%     G12  = Gm·(1+η·Vf)/(1−η·Vf)                [Halpin-Tsai, ξ=1]
%     nu12 = nuf·Vf + num·(1−Vf)                  [ROM]
%
%   GA90R woven corrections (kc=0.92, Naik & Shembekar 1992):
%     E1_GA = E2_GA = 0.5·(E1+E2)·kc
%     nu12_GA = 2·nu12·E2 / (E1+E2)
%     G12_GA  = G12·kc
%
%   Reference: Weisshaar (1981) aeroelastic tailoring; Jones (1999) CLT §2.84.

%% ── Default inputs ─────────────────────────────────────────────────────────
if nargin < 1 || isempty(beta_deg), beta_deg = 20;   end
if nargin < 2 || isempty(Vf),       Vf       = 0.50; end
if nargin < 3,                       matprops = [];   end

%% ── Material constants ─────────────────────────────────────────────────────
% T700 carbon fibre
Ef1   = 230e9;   Ef2  = 15e9;   Gf12 = 27e9;   nuf = 0.20;
rho_f = 1800.0;                                          % kg/m³

% Cured epoxy (LY1564 / Huntsman)
Em    = 3.5e9;   num  = 0.35;   Gm   = Em / (2*(1+num));
rho_m_mat = 1200.0;                                      % kg/m³

% Fabric specifications
FAW_B = 0.300;   % kg/m²  CARBONODB300 ±45 NCF biaxial
FAW_G = 0.302;   % kg/m²  CARBONOGA90R woven 0/90
kc    = 0.92;    % crimp knockdown for woven (Naik & Shembekar 1992)

%% ── Halpin-Tsai micromechanics ─────────────────────────────────────────────
E1   = Ef1*Vf + Em*(1-Vf);               % ROM, fibre-dominated [Pa]
nu12 = nuf*Vf  + num*(1-Vf);             % ROM
E2   = htProp(Ef2,  Em, Vf, 2.0);        % Halpin-Tsai, ξ=2
G12  = htProp(Gf12, Gm, Vf, 1.0);        % Halpin-Tsai, ξ=1

% Apply caller overrides (used by sensitivity study)
if ~isempty(matprops)
    if isfield(matprops, 'E1'),   E1   = matprops.E1;   end
    if isfield(matprops, 'E2'),   E2   = matprops.E2;   end
    if isfield(matprops, 'G12'),  G12  = matprops.G12;  end
    if isfield(matprops, 'nu12'), nu12 = matprops.nu12; end
end

%% ── GA90R woven corrections ────────────────────────────────────────────────
E1_GA   = 0.5*(E1 + E2)*kc;
E2_GA   = E1_GA;
nu12_GA = 2*nu12*E2 / (E1 + E2);
G12_GA  = G12*kc;

%% ── Ply thicknesses from FAW and Vf ────────────────────────────────────────
% FAW is split equally between the +45 and -45 sublayers of each DB300 ply.
t_DB = FAW_B / 2 / (rho_f * Vf);   % [m] single UD sublayer
t_GA = FAW_G     / (rho_f * Vf);   % [m] GA90R woven ply

%% ── AR1 half-stack sequence ────────────────────────────────────────────────
%   1 = DB300 biaxial (expands to 2 sublayers: +45+β, -45+β at t_DB each)
%   0 = GA90R woven   (expands to 1 sublayer:  β      at t_GA)
half_type = [1, 0, 0, 1, 0, 1, 0, 0, 1];  % 4 DB300 + 5 GA90R per half-stack

%% ── Build full symmetric sublayer table ────────────────────────────────────
%   Columns: [angle_deg, thickness_m, is_DB]
subs      = buildSublayers(half_type, beta_deg, t_DB, t_GA);
subs_full = [subs; flipud(subs)];   % mirror bottom half → symmetric laminate

t_total = sum(subs_full(:, 2));     % total thickness [m]
n_sub   = size(subs_full, 1);

%% ── CLT D-matrix integration  (z from −t/2 to +t/2) ───────────────────────
D = zeros(3, 3);
z = -t_total / 2;
for k = 1:n_sub
    ang_k = subs_full(k, 1);
    t_k   = subs_full(k, 2);
    if subs_full(k, 3)   % DB300 — UD Qbar
        Qb = computeQbar(E1, E2, G12, nu12, ang_k);
    else                  % GA90R — woven Qbar
        Qb = computeQbar(E1_GA, E2_GA, G12_GA, nu12_GA, ang_k);
    end
    z0 = z;   z1 = z + t_k;
    D  = D + Qb * (z1^3 - z0^3) / 3;
    z  = z1;
end

D_flex = D;   % [3×3] bending stiffness [N·m]

%% ── Laminate density ────────────────────────────────────────────────────────
rho_lam = rho_f*Vf + rho_m_mat*(1-Vf);   % rule of mixtures [kg/m³]

%% ── D_info struct ───────────────────────────────────────────────────────────
D_info.D11_Nm    = D(1,1);
D_info.D22_Nm    = D(2,2);
D_info.D12_Nm    = D(1,2);
D_info.D66_Nm    = D(3,3);
D_info.D16_Nm    = D(1,3);
D_info.D26_Nm    = D(2,3);
D_info.t_mm      = t_total * 1e3;
D_info.rho_lam   = rho_lam;
D_info.E1_GPa    = E1  / 1e9;
D_info.E2_GPa    = E2  / 1e9;
D_info.G12_GPa   = G12 / 1e9;
D_info.nu12      = nu12;
D_info.t_DB_mm   = t_DB * 1e3;
D_info.t_GA_mm   = t_GA * 1e3;
D_info.beta_deg  = beta_deg;
D_info.Vf        = Vf;
end


%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function y = htProp(Ep, Em_, Vf, xi)
% Halpin-Tsai interpolation: E = Em*(1 + xi*eta*Vf) / (1 - eta*Vf)
%   eta = (Ep/Em - 1) / (Ep/Em + xi)
eta = (Ep/Em_ - 1) / (Ep/Em_ + xi);
y   = Em_ * (1 + xi*eta*Vf) / (1 - eta*Vf);
end


function subs = buildSublayers(half_type, beta_deg, t_DB, t_GA)
% Expand the half-stack type sequence into a sublayer list.
%   DB300 biaxial → 2 sublayers: [+45+β, -45+β], each of thickness t_DB
%   GA90R woven   → 1 sublayer:  [β], thickness t_GA
% Returns [n_sub × 3]: [angle_deg, thickness_m, is_DB]
n    = numel(half_type);
subs = zeros(0, 3);
for k = 1:n
    if half_type(k) == 1   % DB300 biaxial pair
        subs(end+1, :) = [+45 + beta_deg, t_DB, 1]; %#ok<AGROW>
        subs(end+1, :) = [-45 + beta_deg, t_DB, 1]; %#ok<AGROW>
    else                   % GA90R woven
        subs(end+1, :) = [     beta_deg,  t_GA, 0]; %#ok<AGROW>
    end
end
end


function Qb = computeQbar(E1, E2, G12, nu12, theta_deg)
% Transformed reduced stiffness matrix [3×3] in Voigt notation.
% Qbar(θ) maps [ε₁₁, ε₂₂, 2ε₁₂] → [σ₁₁, σ₂₂, σ₁₂].
% Jones (1999) §2.84.
nu21  = nu12 * E2 / E1;
denom = 1 - nu12*nu21;
Q11 = E1  / denom;
Q22 = E2  / denom;
Q12 = nu12*E2 / denom;
Q66 = G12;

c  = cosd(theta_deg);   s  = sind(theta_deg);
c2 = c^2;               s2 = s^2;   cs = c*s;

Qb11 = Q11*c2^2 + 2*(Q12+2*Q66)*s2*c2 + Q22*s2^2;
Qb22 = Q11*s2^2 + 2*(Q12+2*Q66)*s2*c2 + Q22*c2^2;
Qb12 = (Q11+Q22-4*Q66)*s2*c2 + Q12*(c2^2+s2^2);
Qb66 = (Q11+Q22-2*Q12-2*Q66)*s2*c2 + Q66*(s2^2+c2^2);
Qb16 = (Q11-Q12-2*Q66)*c2*cs - (Q22-Q12-2*Q66)*s2*cs;
Qb26 = (Q11-Q12-2*Q66)*s2*cs - (Q22-Q12-2*Q66)*c2*cs;

Qb = [Qb11, Qb12, Qb16;
      Qb12, Qb22, Qb26;
      Qb16, Qb26, Qb66];
end
