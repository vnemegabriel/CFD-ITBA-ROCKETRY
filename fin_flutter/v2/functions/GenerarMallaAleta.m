function mesh = GenerarMallaAleta(cr, ct, span, sweep_rad, nx, ny, varargin)
% GenerarMallaAleta  Structured Q4 mesh for a swept trapezoidal fin.
%
%   mesh = GenerarMallaAleta(cr, ct, span, sweep_rad, nx, ny)
%   mesh = GenerarMallaAleta(..., 'plot', true)
%
%   cr        : root chord [m]
%   ct        : tip chord  [m]
%   span      : semi-span  [m]  (fin height, perpendicular to body)
%   sweep_rad : leading-edge sweep angle [rad]
%   nx        : number of chordwise element divisions
%   ny        : number of spanwise  element divisions
%
%   mesh.nodes   : [(nx+1)*(ny+1) × 3] node XYZ coordinates [m]
%   mesh.connect : [nx*ny × 4] Q4 element connectivity (1-based)
%
%   Coordinate convention:
%     X : chordwise  (LE → TE)
%     Y : spanwise   (root=0 → tip=span)
%     Z : normal out of fin plane (= 0 for flat fin)
%
%   Node numbering: column-major over (xi, eta) grid, xi chordwise, eta spanwise.
%   Element nodes ordered counter-clockwise: [n1 n2 n3 n4] where
%     n1=bottom-left, n2=bottom-right, n3=top-right, n4=top-left.

p = inputParser;
addOptional(p, 'plot', false);
parse(p, varargin{:});
doPlot = p.Results.plot;

xi_v  = linspace(0, 1, nx + 1);    % chordwise parameter [0,1]
eta_v = linspace(0, 1, ny + 1);    % spanwise  parameter [0,1]

nNodes = (nx + 1) * (ny + 1);
nodes  = zeros(nNodes, 3);

k = 0;
for j = 1 : length(eta_v)          % loop spanwise (eta = y/span)
    y   = eta_v(j) * span;
    c_y = cr + (ct - cr) * eta_v(j);          % local chord at this span station
    x0  = y * tan(sweep_rad);                  % LE x-offset due to sweep
    for i = 1 : length(xi_v)       % loop chordwise
        k = k + 1;
        nodes(k, :) = [x0 + xi_v(i) * c_y,  y,  0];
    end
end

% Connectivity: Q4 elements, counter-clockwise node ordering
nEle    = nx * ny;
connect = zeros(nEle, 4);
e = 0;
for j = 1 : ny
    for i = 1 : nx
        e  = e + 1;
        n1 = (j - 1) * (nx + 1) + i;
        n2 = n1 + 1;
        n3 = n2 + (nx + 1);
        n4 = n1 + (nx + 1);
        connect(e, :) = [n1, n2, n3, n4];
    end
end

mesh.nodes   = nodes;
mesh.connect = connect;

if doPlot
    figure('Color','w','Name','Fin Q4 mesh');
    patch('Faces', connect, 'Vertices', nodes, ...
          'FaceColor', [0.85 0.92 1], 'EdgeColor', 'k', 'LineWidth', 0.8);
    axis equal; view(2); grid on;
    xlabel('X – chordwise [m]'); ylabel('Y – spanwise [m]');
    title(sprintf('Swept fin mesh  %d×%d Q4 elements', nx, ny));
end
end
