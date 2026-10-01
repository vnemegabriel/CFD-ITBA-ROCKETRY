function [rho, a, T, P] = isaAtmosphere(h_m)
% isaAtmosphere  International Standard Atmosphere (ISA) model.
%
%   [rho, a, T, P] = isaAtmosphere(h_m)
%
%   h_m  : geometric altitude [m], scalar or vector, 0–32 000 m
%
%   rho  : air density       [kg/m³]
%   a    : speed of sound    [m/s]
%   T    : temperature       [K]
%   P    : pressure          [Pa]
%
%   Layers implemented:
%     0 – 11 000 m  : troposphere       (lapse -6.5 K/km)
%     11 000 – 20 000 m : lower stratosphere (isothermal, T=216.65 K)
%     20 000 – 32 000 m : upper stratosphere (lapse +1.0 K/km)
%
%   Fully vectorised: h_m may be any shape array.

% ---- Constants ----
R     = 287.058;   % [J/(kg·K)] specific gas constant for dry air
gamma = 1.4;       % heat capacity ratio
g0    = 9.80665;   % [m/s²] standard gravity

% ---- Sea-level reference ----
T0  = 288.15;      % [K]
P0  = 101325.0;    % [Pa]
rho0 = 1.225;      % [kg/m³]

% ---- Output pre-allocation ----
sz  = size(h_m);
T   = zeros(sz);
P   = zeros(sz);
rho = zeros(sz);

% ---- Troposphere: 0 – 11 000 m ----
L1    = -6.5e-3;   % lapse rate [K/m]
T11   = T0 + L1 * 11000;         % 216.65 K
P11   = P0 * (T11/T0)^(-g0/(L1*R));

mask1 = h_m <= 11000;
T(mask1)   = T0 + L1 * h_m(mask1);
P(mask1)   = P0 * (T(mask1) / T0) .^ (-g0 / (L1 * R));

% ---- Lower stratosphere: 11 000 – 20 000 m ----
T20   = T11;
P20   = P11 * exp(-g0 * (20000 - 11000) / (R * T11));

mask2 = h_m > 11000 & h_m <= 20000;
T(mask2)   = T11;
P(mask2)   = P11 * exp(-g0 * (h_m(mask2) - 11000) / (R * T11));

% ---- Upper stratosphere: 20 000 – 32 000 m ----
L3    = +1.0e-3;   % lapse rate [K/m]
mask3 = h_m > 20000;
T(mask3)   = T20 + L3 * (h_m(mask3) - 20000);
P(mask3)   = P20 * (T(mask3) / T20) .^ (-g0 / (L3 * R));

% ---- Derived quantities ----
rho = P ./ (R .* T);
a   = sqrt(gamma .* R .* T);
end
