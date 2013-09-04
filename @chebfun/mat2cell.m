function f = mat2cell(f, M, N)
%MAT2CELL   Convert an array-valued CHEBFUN to a cell array of CHEBFUN objects.
%   G = MAT2CELL(F, C) breaks up the array-valued CHEBFUN F into a single row
%   cell array G of CHEBFUN objects. C is the vector of column sizes and must
%   sum to M, the number of columns of F. The elements of C determine the size
%   of each cell in G so that
%               SIZE(C{I},2) == C(I), for I = 1:LENGTH(C)
%
%   G = MAT2CELL(F) assumes C = ones(1, COL).
%
%   G = MAT2CELL(F, M, N) is similar to above, but allows three input arguments
%   so as to be consistent with the built in MAT2CELL function. Here N takes the
%   role of C above, and M should be the scalar value 1.
%
% Example:
%   f = chebfun(@(x) [sin(x), cos(x), exp(x), x], [0, pi])
%   g = mat2cell(f, 1, [1, 2, 1])

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return an empty result:
if ( isempty(f) )
    return
end

% Get the size of the CHEBFUN array:
numFuns = numel(f.funs);
numCols = size(f.funs{1}, 2);

% Parse the inputs:
if ( nargin == 1 )
    M = 1;
    N = ones(1, numCols);
elseif ( nargin == 2 ) 
    N = M;
    M = 1;
end

% Check dimensions:
if ( ~isscalar(M) || M ~= 1 || sum(N) ~= numCols )
    error('CHEBFUN:CHEBFUN:mat2cell:size', ...
        ['Input arguments, M and N, must sum to each dimension of the', ...
        ' input size, [1,%d].'], numCols);
end

% Call MAT2CELL of the FUNs:
cellFuns = cell(numFuns, numel(N));
for k = 1:numFuns
    cellFuns(k,:) = mat2cell(f.funs{k}, M, N);
end

% Create a new CHEBFUN from each column of FUNs;
f = cell(1, numel(N));
for k = 1:numel(N)
    f{k} = chebfun(cellFuns(:,k));
end

end