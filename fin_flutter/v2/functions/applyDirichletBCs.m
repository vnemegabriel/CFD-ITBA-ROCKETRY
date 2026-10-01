function [K_red, M_red, freeDOFs] = applyDirichletBCs(K, M, fixedDOFs)
% applyDirichletBCs  Reduce K and M by eliminating clamped DOFs.
%
%   [K_red, M_red, freeDOFs] = applyDirichletBCs(K, M, fixedDOFs)
%
%   K, M       : full sparse [nDOF × nDOF] matrices
%   fixedDOFs  : vector of DOF indices to clamp (1-based)
%
%   K_red, M_red : reduced sparse matrices (free DOFs only)
%   freeDOFs     : indices of free DOFs in the original system

allDOFs  = (1 : size(K, 1))';
freeDOFs = setdiff(allDOFs, fixedDOFs(:));

K_red = K(freeDOFs, freeDOFs);
M_red = M(freeDOFs, freeDOFs);
end
