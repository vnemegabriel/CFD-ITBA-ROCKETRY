%% mainFlutterSolver.m
%  Supersonic composite fin flutter analysis
%  FEM (Q4 Mindlin shell) + 2nd-order piston theory + p-k flutter solver
%
%  Pipeline:
%    1. Load laminate properties from data/lam.json
%    2. Load flight data from data/flight_data.csv → filter supersonic points
%    3. Build Q4 mesh for swept trapezoidal fin
%    4. Assemble K, M → modal analysis (6 modes)
%    5. Compute GAF at critical Mach via piston theory
%    6. p-k flutter/divergence solve over all supersonic flighxt points
%    7. Report V_flutter and V_div, plot stability envelope
%
%  Run from the supersonic_fin_flutter_matlab/ directory.

clear; clc; close all;
BASE = fileparts(mfilename('fullpath'));
addpath(fullfile(BASE, 'functions'));
fprintf('=== Supersonic Fin Flutter Solver ===\n\n');

%% -----------------------------------------------------------------------
%% 1. Load laminate (data/lam.json)
%% -----------------------------------------------------------------------
lamFile = fullfile(BASE, 'data', 'lam8mm.json'); %Cargar el lam que se quiera analizar
if ~isfile(lamFile)
    error('lam.json not found at %s', lamFile);
end
lam  = jsondecode(fileread(lamFile));
D    = lam.tailored_beta;

D26    = 0;
if isfield(D, 'D26_Nm'), D26 = D.D26_Nm; end
D_flex = [D.D11_Nm, D.D12_Nm, D.D16_Nm;
          D.D12_Nm, D.D22_Nm, D26;
          D.D16_Nm, D26,      D.D66_Nm];

t     = lam.flutter_input.t_mm  * 1e-3;     % shell thickness [m]
rho_m = lam.flutter_input.rho_mat_kgm3;     % material density [kg/m³]

% Isotropic-equivalent for membrane/shear DOFs (derived from D66).
% NOTE: material.E is NOT used by CalcularRigidezQLLL — that function
% re-derives E_eff = 12*D66/t^3 internally from D_flex_3x3.  material.nu
% IS used for the membrane and transverse shear constitutive matrices.
G_eff  = 12 * D.D66_Nm / t^3;
nu_eff = 0.3;
geometry.t  = t;
material.E  = G_eff;   % informational only; overridden inside CalcularRigidezQLLL
material.nu = nu_eff;

fprintf('Laminate (beta=5 tailored, T700/Epoxy AR1):\n');
fprintf('  D11=%.2f  D22=%.2f  D66=%.2f  D16=%.3f  [N·m]\n', ...
        D.D11_Nm, D.D22_Nm, D.D66_Nm, D.D16_Nm);
fprintf('  t=%.2f mm   rho=%.0f kg/m3   E_eff=%.2f GPa\n\n', ...
        t*1e3, rho_m, G_eff/1e9);

%% -----------------------------------------------------------------------
%% 2. Load flight data, filter supersonic points (M >= 1.05)
%% -----------------------------------------------------------------------
csvFile = fullfile(BASE, 'data', 'flight_data.csv');
if ~isfile(csvFile)
    error('flight_data.csv not found at %s', csvFile);
end
opts = detectImportOptions(csvFile, 'CommentStyle', '#');
opts.VariableNames = {'time_s', 'altitude_ft', 'Vz_ms'};
tbl = readtable(csvFile, opts);

h_m = tbl.altitude_ft * 0.3048;
V   = abs(tbl.Vz_ms);

flightPts = struct([]);
for i = 1:height(tbl)
    [rho_i, a_i, ~, ~] = isaAtmosphere(h_m(i));
    M_i = V(i) / a_i;
    if M_i >= 1.05
        fp.Mach  = M_i;   fp.a    = a_i;
        fp.rho   = rho_i; fp.U    = V(i);
        fp.q_inf = 0.5 * rho_i * V(i)^2;
        fp.h_m   = h_m(i); fp.time = tbl.time_s(i);
        if isempty(flightPts), flightPts = fp;
        else, flightPts(end+1) = fp; 
        end
    end
end
if isempty(flightPts)
    error('No supersonic points found in flight_data.csv.');
end

[~, i_crit] = max([flightPts.q_inf]);
fp_crit = flightPts(i_crit);

fprintf('Flight data: %d supersonic points  Mach %.2f-%.2f  h %.0f-%.0f m\n', ...
        numel(flightPts), min([flightPts.Mach]), max([flightPts.Mach]), ...
        min([flightPts.h_m]), max([flightPts.h_m]));
fprintf('Critical point: Mach=%.3f  q=%.0f Pa  h=%.0f m  t=%.1f s\n\n', ...
        fp_crit.Mach, fp_crit.q_inf, fp_crit.h_m, fp_crit.time);



%% -----------------------------------------------------------------------
%% 3. Generate Q4 mesh
%% -----------------------------------------------------------------------
cr        = 0.300;   ct   = 0.150;    % root/tip chord [m]
span      = 0.160;   sweep_deg = 57.4;  % span [m], LE sweep [deg]
nx        = 30;      ny   = 16;

mesh      = GenerarMallaAleta(cr, ct, span, deg2rad(sweep_deg), nx, ny);

%% 3.1 BCs
saveFigs = true;          % save figures to figDir when true

% Fully fixed base
% rootNodes = find(mesh.nodes(:, 2) < 1e-9);
% fixedDOFs = reshape((rootNodes - 1) * 6 + (1:6), 1, []);
% bc_tag    = 'full';

% First two and last two root nodes clamped (LE and TE corners only)
rootNodes   = find(mesh.nodes(:, 2) < 1e-9);
cornerNodes = [rootNodes(1:4); rootNodes(end-3:end)];
fixedDOFs   = reshape((cornerNodes - 1) * 6 + (1:6), 1, []);
bc_tag      = 'corners';

% Middle root nodes clamped – X in [120 mm, 180 mm], 6 chordwise elements
% rootNodes = find(mesh.nodes(:, 2) < 1e-9);
% midNodes  = rootNodes(mesh.nodes(rootNodes, 1) >= 0.120 & ...
%                       mesh.nodes(rootNodes, 1) <= 0.180);
% fixedDOFs = reshape((midNodes - 1) * 6 + (1:6), 1, []);
% bc_tag    = 'mid120_180';

% LE/TE corners + mid patch X in [150 mm, 180 mm]
% rootNodes   = find(mesh.nodes(:, 2) < 1e-9);
% cornerNodes = [rootNodes(1:4); rootNodes(end-3:end)];
% midNodes    = rootNodes(mesh.nodes(rootNodes, 1) >= 0.150 & ...
%                         mesh.nodes(rootNodes, 1) <= 0.180);
% fixedDOFs   = reshape((unique([cornerNodes; midNodes]) - 1) * 6 + (1:6), 1, []);
% bc_tag      = 'corners_mid150_180';


nFixedNodes = length(fixedDOFs) / 6;
divL = repmat('-', 1, 58);
fprintf('\n%s\n', divL);
fprintf('  BC   : %s\n', bc_tag);
fprintf('  Mesh : %d nodes  |  %d x %d Q4 elements\n', ...
        size(mesh.nodes,1), nx, ny);
fprintf('  BCs  : %d root nodes clamped  (%d DOFs fixed)\n', ...
        nFixedNodes, length(fixedDOFs));
fprintf('%s\n\n', divL);

%% -----------------------------------------------------------------------
%% 4. FEM assembly and modal analysis
%% -----------------------------------------------------------------------
K = assembleGlobalStiffness(mesh, geometry, material, D_flex);
M = assembleGlobalMass(mesh, rho_m, t);


[K_red, M_red, freeDOFs] = applyDirichletBCs(K, M, fixedDOFs);

nModes = 6;
[Phi_red, omega_n] = modalAnalysis(K_red, M_red, nModes);
f_n = omega_n / (2*pi);

Phi_full = zeros(size(K,1), nModes);
Phi_full(freeDOFs, :) = Phi_red;

fprintf('  Frequencies [Hz] :'); fprintf(' %8.2f', f_n); fprintf('\n\n');


%% -----------------------------------------------------------------------
%% 5. GAF at critical flight condition
%% -----------------------------------------------------------------------
% Reference semi-chord (MAC/2) for reduced-frequency k = omega*b_ref/U
b_ref  = (cr^2 + cr*ct + ct^2) / (3*(cr + ct));
k_vals = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0];

fprintf('Computing GAF (Mach=%.3f)... ', fp_crit.Mach);
Q_k      = pistonTheoryGAF(mesh, Phi_full, fp_crit.Mach, fp_crit.q_inf, ...
                            fp_crit.a, k_vals, sweep_deg);
Q_k_norm = Q_k / fp_crit.q_inf;   % per unit dynamic pressure [1/Pa]
fprintf('done\n\n');

%% -----------------------------------------------------------------------
%% 6. p-k flutter/divergence solve
%% -----------------------------------------------------------------------
fprintf('Running p-k solver over %d flight points...\n', numel(flightPts));
[V_fl, V_div, pkRes] = pkSolveFlutter(omega_n, Q_k_norm, k_vals, b_ref, flightPts);
fprintf('done\n\n');

% ── Summary ──────────────────────────────────────────────────────────────
V_flight = [flightPts.U]';
mach_vec = [flightPts.Mach]';

fl_margin  = V_fl  ./ V_flight;
div_margin = V_div ./ V_flight;

[min_fl_margin,  i_fl]  = min(fl_margin);
[min_div_margin, i_div] = min(div_margin);

% ── Structured results block (copy-paste ready for comparison) ────────────
divE = repmat('=', 1, 58);
fprintf('%s\n', divE);
fprintf('  RESULTS  —  BC: %s\n', bc_tag);
fprintf('%s\n', divE);
fprintf('  Fixed nodes      : %d   Fixed DOFs: %d\n', nFixedNodes, length(fixedDOFs));
fprintf('  Frequencies [Hz] :'); fprintf(' %8.2f', f_n); fprintf('\n');
fprintf('  Flutter  (p-k)   : ');
if all(isinf(V_fl))
    fprintf('STABLE  (all %d pts, g_struct=%.2f)\n', numel(flightPts), pkRes.g_struct);
else
    fprintf('V_fl  = %6.0f m/s   Mach %.3f   margin %.2fx\n', ...
            V_fl(i_fl), flightPts(i_fl).Mach, min_fl_margin);
end
fprintf('  Divergence (k=0) : ');
if all(isinf(V_div))
    fprintf('STABLE  (all %d pts)\n', numel(flightPts));
else
    fprintf('V_div = %6.0f m/s   Mach %.3f   margin %.2fx\n', ...
            V_div(i_div), flightPts(i_div).Mach, min_div_margin);
end
fprintf('%s\n\n', divE);

% ── Append one row to bc_comparison.csv (accumulates across BC runs) ──────
csvOut     = fullfile(BASE, 'bc_comparison.csv');
fl_stable  = double(all(isinf(V_fl)));
div_stable = double(all(isinf(V_div)));
V_fl_out   = min(V_fl(isfinite(V_fl)));   if isempty(V_fl_out),  V_fl_out  = -1; end
V_div_out  = min(V_div(isfinite(V_div))); if isempty(V_div_out), V_div_out = -1; end
fl_mar_out  = min(fl_margin(isfinite(fl_margin)));   if isempty(fl_mar_out),  fl_mar_out  = -1; end
div_mar_out = min(div_margin(isfinite(div_margin))); if isempty(div_mar_out), div_mar_out = -1; end
if ~isfile(csvOut)
    fid = fopen(csvOut, 'w');
    fprintf(fid, 'bc_tag,fixed_nodes,fixed_DOFs,f1_Hz,f2_Hz,f3_Hz,f4_Hz,f5_Hz,f6_Hz,flutter_stable,V_fl_min_ms,fl_margin,div_stable,V_div_min_ms,div_margin\n');
    fclose(fid);
end
fid = fopen(csvOut, 'a');
fprintf(fid, '%s,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%.1f,%.3f,%d,%.1f,%.3f\n', ...
    bc_tag, nFixedNodes, length(fixedDOFs), ...
    f_n(1), f_n(2), f_n(3), f_n(4), f_n(5), f_n(6), ...
    fl_stable, V_fl_out, fl_mar_out, div_stable, V_div_out, div_mar_out);
fclose(fid);

%% -----------------------------------------------------------------------
%% 7. Figures
%% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
% Results folder naming depending the thickness
% -----------------------------------------------------------------------
th_mm = lam.inputs.target_thickness_mm;          % thickness in mm (may be numeric or string)
if isnumeric(th_mm)
    th_str = num2str(th_mm);
else
    th_str = th_mm;
end
th_str = strrep(th_str, '.', '_');              % replace decimal point with underscore
figDir = fullfile(BASE, ['results' th_str 'mm']);
if ~exist(figDir, 'dir')
    mkdir(figDir);
end
[mach_s, si] = sort(mach_vec);
% -----------------------------------------------------------------------

% ── Fig 1: Fin mesh ───────────────────────────────────────────────────────
figure('Color','w','Visible','on');
patch('Faces', mesh.connect, 'Vertices', mesh.nodes(:,1:2), ...
      'FaceColor',[0.88 0.93 1],'EdgeColor',[0.4 0.4 0.4],'LineWidth',0.5);
axis equal tight; grid on; box on;
xlabel('Chordwise x [m]'); ylabel('Spanwise y [m]');
title(sprintf('Q4 Mindlin shell mesh — %d×%d elements\ncr=%.0f mm, ct=%.0f mm, span=%.0f mm, \\Lambda=%.1f°', ...
              nx, ny, cr*1e3, ct*1e3, span*1e3, sweep_deg));
drawnow;
if saveFigs
    saveas(gcf, fullfile(figDir, [bc_tag '_mesh.png']));
end

% ── Fig 2: Mode shapes (patch, interpolated colour) ───────────────────────
nN    = size(mesh.nodes, 1);
wDOFs = (0:nN-1)' * 6 + 3;       % out-of-plane DOF per node

figure('Color','w','Visible','on','Position',[100 100 900 700]);
for mi = 1:min(4, nModes)
    subplot(2, 2, mi);
    w_vals = Phi_full(wDOFs, mi);
    % Normalise for display
    w_vals = w_vals / max(abs(w_vals) + eps);
    patch('Faces', mesh.connect, 'Vertices', mesh.nodes(:,1:2), ...
          'FaceVertexCData', w_vals, 'FaceColor','interp', ...
          'EdgeColor',[0.3 0.3 0.3], 'EdgeAlpha', 0.25, 'LineWidth', 0.3);
    colormap(gca, 'parula'); 
    caxis([-1 1]);
    colorbar('Location','eastoutside');
    axis equal tight; grid off; box on;
    title(sprintf('Mode %d — f_n = %.1f Hz', mi, f_n(mi)));
    xlabel('x [m]'); ylabel('y [m]');
end
sgtitle('Mass-normalised mode shapes (out-of-plane w, normalised to max=1)');
drawnow;
if saveFigs
    saveas(gcf, fullfile(figDir, [bc_tag '_mode_shapes.png']));
end

% ── Fig 3: Velocity envelope — V_flight vs V_flutter / V_div ──────────────
time_s  = [flightPts.time]';
V_flt_s = V_fl(si);
V_div_s = V_div(si);
V_flt_s = min(V_flt_s, 5000);   % cap Inf for plotting
V_div_s = min(V_div_s, 5000);

figure('Color','w','Visible','on','Position',[100 100 900 500]);
hold on;
% Safe corridor fill (between flight speed and minimum critical speed)
V_crit = min(V_flt_s, V_div_s);
fill([mach_s; flipud(mach_s)], [V_flight(si); flipud(V_crit)], ...
     [0.85 1 0.85], 'EdgeColor','none', 'FaceAlpha', 0.4, ...
     'DisplayName','Safe corridor');

plot(mach_s, V_flight(si),    'k-',  'LineWidth', 2.0, 'DisplayName','Flight speed');
plot(mach_s, V_flt_s,         'b--', 'LineWidth', 1.8, 'DisplayName','Flutter speed (p-k)');
plot(mach_s, V_div_s,         'r--', 'LineWidth', 1.8, 'DisplayName','Divergence speed');

% Annotate if capped
if any(isinf(V_fl))
    text(0.98, 0.95, 'V_{flutter} = \infty  (wash-out design)', ...
         'Units','normalized','HorizontalAlignment','right', ...
         'Color','b','FontSize',9,'FontWeight','bold');
end
if any(isinf(V_div))
    text(0.98, 0.88, 'V_{div} = \infty  (all div. eigenvalues \leq 0)', ...
         'Units','normalized','HorizontalAlignment','right', ...
         'Color','r','FontSize',9,'FontWeight','bold');
end

hold off; grid on; box on;
xlabel('Mach number  (M \geq 1.05, piston theory valid)');
ylabel('Velocity [m/s]');
title(sprintf(['Aeroelastic velocity envelope  —  \\Lambda=%.1f° swept fin + D_{16}=%.2f N·m\n' ...
               '%d/%d supersonic points flutter-safe   %d/%d divergence-safe'], ...
              sweep_deg, D_flex(1,3), ...
              sum(~isinf(V_fl) & V_fl > V_flight | isinf(V_fl)), numel(flightPts), ...
              sum(~isinf(V_div) & V_div > V_flight | isinf(V_div)), numel(flightPts)));
legend('Location','best','FontSize',9);
% drawnow;
if saveFigs
    saveas(gcf, fullfile(figDir, [bc_tag '_velocity_envelope.png']));
end

% ── Fig 4: p-k damping history at critical flight point ───────────────────
gam_crit = pkRes.gam_hist{i_crit};
q_crit   = pkRes.q_hist{i_crit};

figure('Color','w','Visible','on','Position',[100 100 900 450]);
cmap = lines(nModes);
hold on;
for mi = 1:nModes
    plot(q_crit/1e3, gam_crit(mi,:), '-', 'LineWidth', 1.5, ...
         'Color', cmap(mi,:), 'DisplayName', sprintf('Mode %d (%.0f Hz)', mi, f_n(mi)));
end
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', '\gamma = 0  (flutter onset)');
xline(fp_crit.q_inf/1e3, 'k:', 'LineWidth', 1.2, 'DisplayName', 'q_{flight}');
hold off; grid on; box on;
xlabel('Dynamic pressure q  [kPa]');
ylabel('Growth rate \gamma = Re(p)  [rad/s]');
title(sprintf('p-k damping history — Mach=%.3f, h=%.0f m, g_{struct}=%.2f', ...
              fp_crit.Mach, fp_crit.h_m, pkRes.g_struct));
legend('Location','best','FontSize',8);
% drawnow;
if saveFigs
    saveas(gcf, fullfile(figDir, [bc_tag 'pk_damping.png']));
end

%% -----------------------------------------------------------------------
%% 8. Save results
%% -----------------------------------------------------------------------
save(fullfile(figDir, 'flutter.mat'), ...
    'V_fl', 'V_div', 'fl_margin', 'div_margin', ...
    'omega_n', 'f_n', 'Phi_full', 'mesh', ...
    'flightPts', 'fp_crit', 'Q_k', 'Q_k_norm', 'k_vals', 'b_ref', ...
    'pkRes', 'D_flex', 't', 'rho_m');

fprintf('flutter.mat  → %s\n', figDir);
fprintf('Figures      → %s/[figure]_%s.png\n', figDir, bc_tag);
fprintf('Comparison   → %s\n', csvOut);
fprintf('\n=== Done — BC: %s ===\n\n', bc_tag);
