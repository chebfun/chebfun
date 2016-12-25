function coeffs3D = vals2coeffs(F)
%VALS2COEFFS   Convert tensor of values to tensor of coefficients.
%
% See also CHEBFUN3/COEFFS2VALS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[m, n, p] = size(F);

%% Step 1:
% Mode-1 unfolding of F to get a matrix of size m x n*p:
F1 = chebfun3.unfold(F, 1); 

% Apply 1D "vals2coeffs" in the X direction i.e., the 1st dimension of F1:
F1 = chebtech2.vals2coeffs(F1);

% Tensorize F1 back to its original m x n x p size:
F1 = chebfun3.fold(F1, [m, n, p], 1, [2 3]); 
% This is simply F1 = reshape(V1, m, n, p);

%% Step 2:
% Mode-2 unfolding of (the tensorized) F1 to get a matrix of size m x n*p:
F2 = chebfun3.unfold(F1, 2); 

% Apply 1D "vals2coeffs" in the Y direction, i.e. to the 1st direction of
% F2:
F2 = chebtech2.vals2coeffs(F2);

% Tensorize F2 back to its original m x n x p size:
F2 = chebfun3.fold(F2, [m, n, p], 2, [1 3]);

%% Step 3:
% Mode-3 unfolding of (the tensorized) F2 to get a matrix of size p x m*n:
F3 = chebfun3.unfold(F2, 3);

% Now, vals2coeffs is applied in the Z direction.
F3 = chebtech2.vals2coeffs(F3);

% Reshape F3 back to the original m x n x p size:
coeffs3D = chebfun3.fold(F3, [m, n, p], 3, [1 2]);
end