function T = outerProd(A, B)
%OUTERPROD   Outer product 
%   OUTERPROD(A, B) computes outer product of the matrix A and the vector B
%   or a vector A and a matrix B. Both A and B are discrete.
%
%   Example:
%   A = rand(4,1); 
%   B = rand(3,2);
%   T = outerProd(A, B);       % T will be 4 x 1 x 3 x 2.
%
% See also CHEBFUN3/TXM.

% The structure of this code is similar to `ttt.m` from the HTUCKER toolbox
% of Tobler and Kressner.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Determine number of dimensions in the inputs:
ndimsA = ndims(A);
ndimsB = ndims(B);

% Determine sizes of the inputs:
sizeA = size(A);
sizeA(end+1:ndimsA) = 1;
sizeB = size(B);
sizeB(end+1:ndimsB) = 1;

% Vectorize both input matrices:
% Call chebfun3.unfold:
vecA = A(:);
% The same as vecA = chebfun3.unfold(A, [1:ndimsA], []);

% Call chebfun3.unfold:
vecB = B(:);
% The same as vecB = chebfun3.unfold(B, [1:ndimsB], []);

% Calculate matricization of T_mat, i.e., T_mat = vecA * vecB.'
% The outer product is performed in the vector level.
T_mat = vecA * vecB.'; 

% Determine size of the output tensor T:
size_T = [sizeA, sizeB];

% Enforce at least two entries in size_T by possibly adding singleton
% dimensions
size_T(end+1:2) = 1;

% Fold the matrix zMat to the tensor T_mat
T = chebfun3.fold(T_mat, size_T, 1:2, 3:numel(size_T));

end