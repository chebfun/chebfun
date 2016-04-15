function A = unfold(T, rowDims, colDims)
%UNFOLD     Unfold or flatten or matricize a tensor.
%
%   A = UNFOLD(T, ROWDIMS, COLDIMS) returns an unfolding of T, where 
%   ROWDIMS and COLDIMS are integer vectors denoting the dimensions of T 
%   that are merged into the rows and the columns of the matrix,
%   respectively.
%
%   A = UNFOLD(T, ROWDIMS) chooses COLDIMS automatically to contain all 
%   dimensions not contained in ROWDIMS. The entries of COLDIMS are
%   arranged in ascending order.
%
%   Note that UNFOLD(T, 1:ndims(T)) is the same as to T(:).
%
%   Example
%   T = randn(2, 3, 4, 5);
%   A = unfold(T, [2 4]); % returns a 15x8 matrix
%
%   See also CHEBFUN3/FOLD.

% The structure of this code is similar to `matricize.m` from the HTUCKER 
% toolbox of Tobler and Kressner.

if ( nargin == 3 )
    rowDims = rowDims(rowDims <= ndims(T));
    colDims = colDims(colDims <= ndims(T));  
elseif ( nargin == 2 )
    % Put remaining modes onto columns.
    colDims = setdiff(1:ndims(T), rowDims);
else
    error('CHEBFUN:CHEBFUN3:unfold', ...
        'COL_DIMS must be a positive integer vector without duplicate entries.')
end

% Size of tensor T
sizeT = size(T);
if ( ndims(T) < 3 )
    % To work properly also with scalars, vectors and matrices.
    sizeT(3) = 1;
end

% Permute dimensions of T
T = permute(T, [rowDims, colDims]);

% Flatten the tensor T
A = reshape(T, prod(sizeT(rowDims)), prod(sizeT(colDims)) );

end