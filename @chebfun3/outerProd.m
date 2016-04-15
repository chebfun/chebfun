function T = outerProd(A, B)
%OUTERPROD    Outer product of a matrix and a vector, or a vector and a 
%             matrix, all discrete data.
%
%   Example
%   A = rand(4,1); 
%   B = rand(3,2);
%   T = outerProd(A, B);       % T will be 4 x 1 x 3 x 2.
%   See also CHEBFUN3/TXM.

% The structure of this code is similar to `ttt.m` from the HTUCKER toolbox
% of Tobler and Kressner.

% Determine number of dimensions in the inputs:
ndimsA = ndims(A);
ndimsB = ndims(B);

% Determine sizes of the inputs:
sizeA = size(A);
sizeA(end+1:ndimsA) = 1;
sizeB = size(B);
sizeB(end+1:ndimsB) = 1;

% Vectorize both input matrices:
%vecA = chebfun3.unfold(A, [1:ndimsA], []);
vecA = A(:);

%vecB = chebfun3.unfold(B, [1:ndimsB], []);
vecB = B(:);

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