function [K] = assembleGlobalStiffness(mesh, geometry, material, D_flex_3x3)
% assembleGlobalStiffness  Sparse global stiffness matrix for Q4 Mindlin shell.
%
%   K = assembleGlobalStiffness(mesh, geometry, material, D_flex_3x3)
%
%   mesh         : struct with .nodes [nNodes×3] and .connect [nEle×4]
%   geometry.t   : shell thickness [m]
%   material.nu  : effective Poisson ratio (0.3 recommended for composite)
%   D_flex_3x3   : [3×3] CLT anisotropic bending stiffness matrix [N·m]
%
%   K            : sparse [6*nNodes × 6*nNodes] global stiffness matrix

nNodes = size(mesh.nodes, 1);
nDOF   = 6 * nNodes;
K      = sparse(nDOF, nDOF);

nEle   = size(mesh.connect, 1);
nDOF_e = 24;                           % 4 nodes × 6 DOFs

% Pre-allocate triplet storage for sparse assembly
I = zeros(nEle * nDOF_e^2, 1);
J = zeros(nEle * nDOF_e^2, 1);
V = zeros(nEle * nDOF_e^2, 1);
cnt = 0;

for ie = 1:nEle
    nd      = mesh.connect(ie, :);     % [1×4] node indices
    nodes3D = mesh.nodes(nd, :);       % [4×3]

    Ke = CalcularRigidezQLLL(nodes3D, geometry, material, D_flex_3x3, 'selective');

    % Global DOF indices for this element: node-major ordering
    % CalcularRigidezQLLL uses node-major local DOFs:
    %   [u1,v1,w1,θx1,θy1,θz1, u2,..., u4,...,θz4]
    % Transpose the 4×6 index matrix so reshape reads rows (nodes) first.
    dofs = reshape(((nd - 1)' * 6 + (1:6))', 1, []);  % [1×24] node-major

    [gJ, gI] = meshgrid(dofs, dofs);
    r = cnt + 1 : cnt + nDOF_e^2;
    I(r) = gI(:);
    J(r) = gJ(:);
    V(r) = Ke(:);
    cnt  = cnt + nDOF_e^2;
end

K = sparse(I, J, V, nDOF, nDOF);
end
