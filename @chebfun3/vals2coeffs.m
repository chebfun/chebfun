function coeffs3D = vals2coeffs(vals3D)
%VALS2COEFFS     Convert tensor of values to tensor of coefficients
%
%   See also chebfun3/coeffs2vals.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m, n, p] = size(vals3D);

%% Step 1:
% Mode-1 unfolding of F to get a matrix of size m x n*p:
vals1 = chebfun3.unfold(vals3D, 1); 

% Apply 1D "vals2coeffs" in the X direction i.e., the 1st dimension of vals1:
vals1 = chebtech2.vals2coeffs(vals1); 

% Tensorize vals1 back to its original m x n x p size:
vals1 = chebfun3.fold(vals1, [m, n, p], 1, [2 3]); 
% This is simply vals1 = reshape(V1, m, n, p);

%% Step 2:
% Mode-2 unfolding of (the tensorized) vals1 to get a matrix of size m x n*p:
vals2 = chebfun3.unfold(vals1, 2); 

% Apply 1D "vals2coeffs" in the Y direction, i.e. to the 1st direction of
% vals2:
vals2 = chebtech2.vals2coeffs(vals2);

% Tensorize vals2 back to its original m x n x p size:
vals2 = chebfun3.fold(vals2, [m, n, p], 2, [1 3]);

%% Step 3:
% Mode-3 unfolding of (the tensorized) vals2 to get a matrix of size p x m*n:
vals3 = chebfun3.unfold(vals2, 3);

% Now, vals2coeffs is applied in the Z direction.
vals3 = chebtech2.vals2coeffs(vals3);

% Reshape vals3 back to the original m x n x p size:
coeffs3D = chebfun3.fold(vals3, [m, n, p], 3, [1 2]);
end