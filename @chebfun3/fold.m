function T = fold(A, sizeT, rowDims, colDims)
%FOLD   Reshapes a matrix to get a discrete tensor.
%   T = FOLD(A, SIZET, ROWDIMS, COLDIMS) returns the tensor T of size
%   SIZET. This is the reverse of A = UNFOLD(T, ROWDIMS, COLDIMS),
%   where ROWDIMS and COLDIMS contain the dimensions of T that correspond
%   to the rows and to the columns of A, respectively.
%
%   Example 1:
%   A = randn(15, 8);
%   T = chebfun3.fold(A, [2 3 4 5], [2 4]) % returns a 2x3x4x5 tensor
%   
%   Example 2: 
%   T = rand(m,n,p);
%   T2 = chebfun3.unfold(T, 2); % Mode-2 unfolding of T
%
%   Here is the reverse operation:
%   T_new = chebfun3.fold(T2, [m, n, p], 2, [1 3]);
%
% See also CHEBFUN3/UNFOLD.

% The structure of this code is similar to `dematricize.m` from the HTUCKER 
% toolbox of Tobler and Kressner.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

sizeT(end+1:2) = 1;
dim = numel(sizeT);
rowDims = rowDims(rowDims <= dim);
colDims = colDims(colDims <= dim);

% Reshape to tensor
T = reshape(A, sizeT([rowDims, colDims]) );

% Inverse permute
T = ipermute(T, [rowDims, colDims]);

end