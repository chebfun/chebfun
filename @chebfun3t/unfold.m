function A = unfold(T, rowDims, colDims)
%UNFOLD   Unfold or matricize a discrete tensor.
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
%   Example:
%   T = randn(2, 3, 4, 5);
%   A = unfold(T, [2 4]); % returns a 15x8 matrix
%
% See also CHEBFUN3T/FOLD.

% The structure of this code is similar to `matricize.m` from the HTUCKER 
% toolbox of Tobler and Kressner.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 3 )
    % No argument checking, remove trailing singleton dimensions
    rowDims = rowDims(rowDims <= ndims(T));
    colDims = colDims(colDims <= ndims(T));  
elseif ( nargin == 2 )
    colDims = setdiff(1:ndims(T), rowDims);
else
    error('CHEBFUN:CHEBFUN3:unfold', ...
        'COL_DIMS must be a positive integer vector without duplicate entries.')
end

% Save size of original tensor
sizeT = size(T);

% Permute dimensions of tensor x
T = permute(T, [rowDims, colDims]);

% Reshape to matrix
A = reshape(T, prod(sizeT(rowDims)), prod(sizeT(colDims)) );

end