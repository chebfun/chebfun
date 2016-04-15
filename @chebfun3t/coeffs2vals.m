function V = coeffs2vals(C)
%VAL2COEFFS    Convert tensor of Chebyshev coefficients to values.
%
%   V = COEFFS2VALS(C) converts a tensor C of trivariate Chebyshev 
%   coefficients to a tensor of samples V corresponding to values of 
%   sum_i sum_j sum_k C(i,j,k) T_{i-1}(x) T_{j-1}(y) T_{k-1}(z) 
%   at a tensor Chebyshev grid of the same size as C.
%
%   See also chebfun3t/vals2coeffs.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m, n, p] = size(C);

%% Step 1:
% Mode-1 unfolding of F to get a matrix of size m x n*p:
C1 = chebfun3t.unfold(C, 1); 

% Apply 1D "vals2coeffs" in the X direction i.e., the 1st dimension of vals1:
C1 = chebtech2.coeffs2vals(C1); 

% Tensorize vals1 back to its original m x n x p size:
C1 = chebfun3t.fold(C1, [m, n, p], 1, [2 3]); 
% This is simply vals1 = reshape(V1, m, n, p);

%% Step 2:
% Mode-2 unfolding of (the tensorized) vals1 to get a matrix of size m x n*p:
C2 = chebfun3t.unfold(C1, 2); 

% Apply 1D "vals2coeffs" in the Y direction, i.e. to the 1st direction of
% vals2:
C2 = chebtech2.coeffs2vals(C2);

% Tensorize vals2 back to its original m x n x p size:
C2 = chebfun3t.fold(C2, [m, n, p], 2, [1 3]);

%% Step 3:
% Mode-3 unfolding of (the tensorized) vals2 to get a matrix of size p x m*n:
C3 = chebfun3t.unfold(C2, 3);

% Now, vals2coeffs is applied in the Z direction.
C3 = chebtech2.coeffs2vals(C3);

% Reshape vals3 back to the original m x n x p size:
V = chebfun3t.fold(C3, [m, n, p], 3, [1 2]);

end