function [M] = assembleGlobalMass(mesh, rho, t)
% assembleGlobalMass  Sparse lumped mass matrix for Q4 Mindlin shell.
%
%   M = assembleGlobalMass(mesh, rho, t)
%
%   mesh  : struct with .nodes [nNodes×3] and .connect [nEle×4]
%   rho   : material density [kg/m³]
%   t     : shell thickness [m]
%
%   M     : sparse [6*nNodes × 6*nNodes] lumped mass matrix
%           Each nodal translational & rotational DOF gets rho*t*A_elem/4.
%           Rotational inertia: (rho*t^3/12)*A_elem/4 per rotational DOF.

nNodes = size(mesh.nodes, 1);
nDOF   = 6 * nNodes;
nEle   = size(mesh.connect, 1);
nDOF_e = 24;

I = zeros(nEle * nDOF_e, 1);
V = zeros(nEle * nDOF_e, 1);
cnt = 0;

for ie = 1:nEle
    nd   = mesh.connect(ie, :);
    pts  = mesh.nodes(nd, :);          % [4×3]

    % Quadrilateral area: sum of two triangles (exact for general quads)
    t1   = 0.5 * norm(cross(pts(2,:) - pts(1,:), pts(3,:) - pts(1,:)));
    t2   = 0.5 * norm(cross(pts(3,:) - pts(1,:), pts(4,:) - pts(1,:)));
    area = t1 + t2;

    m_trans = rho * t   * area / 4;   % translational mass per node [kg]
    m_rot   = rho * t^3 * area / 48;  % rotational inertia per node [kg·m²]
    %   t^3/12 * area/4 = t^3*area/48

    % Lumped diagonal: [u, v, w, θx, θy, θz] per node
    m_node = [m_trans; m_trans; m_trans; m_rot; m_rot; m_rot];

    dofs = reshape(((nd - 1)' * 6 + (1:6))', [], 1);  % [24×1] node-major
    m_vec = repmat(m_node, 4, 1);                    % [24×1]

    r = cnt + 1 : cnt + nDOF_e;
    I(r) = dofs;
    V(r) = m_vec;
    cnt  = cnt + nDOF_e;
end

M = sparse(I, I, V, nDOF, nDOF);
end
