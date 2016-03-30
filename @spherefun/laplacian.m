function f = laplacian(f) 
%LAPLACIAN   Scalar laplacian of a SPHEREFUN.
%   L = LAPLACIAN(F) computes the spherical Laplacian of the SPHEREFUN F.  
%  
% This is gvien, in intrinsic spherical coordinates as
% 
% laplacian(F) = d_{theta} d_{theta} F + tan(theta) d_{theta} F + 
%                1/sin(theta) + d_{lambda} d_{lambda} F
%
% where 0 <= theta <= pi is the co-latitude variable and 
% -pi <= lambda <= pi  is the longitude variable.
%
% See also GRADIENT, DIVERGENCE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% We implement the computation of the Laplacian using the spherefun/diff
% command.  This function computes tangential derivatives in the x, y, and
% z directions, so that the Laplacian can be implemented as we usually
% think of it in 3D, i.e.
% laplacian(f) = f_xx + f_yy + f_zz.

% Compute the derivatives
fxx = diff(f, 1, 2); 
fyy = diff(f, 2, 2); 
fzz = diff(f, 3, 2);

% Determine which derivative has the most coefficients in its
% representation for the columns and rows.
[mxx,nxx] = length(fxx);
[myy,nyy] = length(fyy);
[mzz,nzz] = length(fzz);

m = max([mxx myy mzz]);
n = max([nxx nyy nzz]);

% Ensure m and n are even:
m = m + mod(m, 2);
n = n + mod(n, 2);

% Construct a spherefun by sampling the derivatives on the same size grid.
f = spherefun(sample(fxx, m, n/2) + sample(fyy, m, n/2) + ...
               sample(fzz, m, n/2));

end
