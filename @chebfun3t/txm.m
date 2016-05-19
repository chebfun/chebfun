function T = txm(T, A, varargin)
%TXM   Modal contraction of a discrete Tensor by a matrix.
%   Y = TXM(T, A, K) computes mode-K multiplication of a tensor T by a 
%   matrix A. Inputs must satisfy the relation 
%                       size(T, K) == size(A, 2)
%   and the output Y satisfies size(Y, K) = size(A, 1).
%
%   Example:
%   T = rand([3, 5, 2, 4]);
%   A = rand(6, 3); 
%   S = txm(T, A, 1);  % Mode-1 contraction of T with A. S will be a 
%                        tensor of size 6 x 5 x 2 x 4.
%
% See also CHEBFUN3T/outerProd.

% The structure of this code is similar to `ttm.m` from the HTUCKER toolbox 
% of Tobler and Kressner.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

contractionMode = varargin{1};
dim = max(ndims(T), max(contractionMode));

% Calculate size (add singleton dimensions if A contains more entries)
size_T = size(T);
size_T(end+1:dim) = 1;

% Unfold the tensor T:
unContractedModes = [1:contractionMode-1, contractionMode+1:dim];  
if ( contractionMode == 1 )
    % T will not be transposed.
    T_mat = chebfun3t.unfold(T, contractionMode, unContractedModes);
else
    % T will be transposed.
    T_mat = chebfun3t.unfold(T, unContractedModes, contractionMode);
end

if( size(A, 2) ~= size_T(contractionMode) )
    error('CHEBFUN:CHEBFUN3:txm', 'Tensor-matrix dimensions do not agree');
end

% Perform the modal contraction:
if ( contractionMode == 1 )    
    T_mat = A * T_mat;
else
    T_mat = T_mat * A.';
end

% Update size of T
if ( contractionMode == 1 )    
    size_T(contractionMode) = size(T_mat, 1);
else
    size_T(contractionMode) = size(T_mat, 2);
end

% Fold T_mat
if ( contractionMode == 1 )    
    T = chebfun3t.fold(T_mat, size_T, contractionMode, ...
        unContractedModes);
else
    T = chebfun3t.fold(T_mat, size_T, unContractedModes, ...
        contractionMode);
end

end