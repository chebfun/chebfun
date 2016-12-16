function g = mat2cell(f, M, N)
%MAT2CELL   Convert an array-valued CLASSICFUN into a cell array of CLASSICFUN objects.
%   G = MAT2CELL(F, C) breaks up the array-valued CLASSICFUN F into a single row cell
%   array G of CLASSICFUN objects. C is the vector of column sizes and must sum to M,
%   the number of columns of F. The elements of C determine the size of each
%   cell in G so that
%               SIZE(C{I},2) == C(I), for I = 1:LENGTH(C)
%
%   G = MAT2CELL(F) assumes C = ones(1, COL).
%
%   G = MAT2CELL(F, M, N) is similar to above, but allows three input arguments
%   so as to be consistent with the built in MAT2CELL function. Here N takes the
%   role of C above, and M should be the scalar value 1.
%
% Example:
%   f = bndfun(@(x) [sin(x), cos(x), exp(x), x], [0, pi])
%   g = mat2cell(f, 1, [1, 2, 1])

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return an empty result:
if ( isempty(f) )
    g = f;
    return
end

% Get the size of the onefun array:
m = size(f.onefun, 2);

% Parse the inputs:
if ( nargin == 1 )
    M = 1;
    N = ones(1, m);
elseif ( nargin == 2 ) 
    N = M;
    M = 1;
end

% Check dimensions:
if ( ~isscalar(M) || M ~= 1)
    error('CHEBFUN:CLASSICFUN:mat2cell:size', ...
        ['Input arguments, M and N, must sum to each dimension of the', ...
        ' input size, [1,%d].'], size(f.onefun, 2));
end

% Call MAT2CELL of the ONEFUN:
oneFuns = mat2cell(f.onefun, M, N);

% Create a cell for storing the new CLASSICFUN objects:
g = cell(1, numel(N));

% Append the data to the new entries in the cell:
for k = 1:numel(N)
    % Create a new CLASSICFUN from the ONEFUN and the domain:
    data.domain = f.domain;
    g{k} = f.make(oneFuns{k}, data);
end

end
