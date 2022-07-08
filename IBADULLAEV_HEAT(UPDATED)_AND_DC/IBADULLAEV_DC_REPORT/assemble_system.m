function [A,b] = assemble_system(type, x, y, k, I, global_idx)
%ASSEMBLE_SYSTEM Summary of this function goes here
%   Detailed explanation goes here
nx = length(x);
ny = length(y);
n = nx*ny;

[C0, C1, C2, C3, C4] = assemble_coefficients(x, y, k, type);
A_ = spdiags([C3(:), C1(:), C0(:), C2(:), C4(:)], ...
              nx+[-nx, -1, 0, 1, nx], n, n+2*nx); % initialize as recangular matrix
A = A_(:, nx+1:end-nx);                         % remove redundant columns
b = assemble_rhs(x, y, global_idx, I);
end

