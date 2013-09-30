function g = mat2cell(f, M, N)
%MAT2CELL   Convert an array-valued CHEBFUN to a cell array of CHEBFUN objects.
%   G = MAT2CELL(F, C) breaks up the array-valued CHEBFUN F into a cell array G
%   of CHEBFUN objects. C is a vector of sizes and must sum to the number of
%   components of F (i.e., the number of columns (rows) of F if F is a column
%   (row) CHEBFUN). The elements of C determine the size of each cell in G so
%   that
%               SIZE(G{I}, 2) == C(I), for I = 1:LENGTH(C)
%   if F is a column CHEBFUN or
%               SIZE(G{I}, 1) == C(I), for I = 1:LENGTH(C)
%   if F is a row CHEBFUN.
%
%   G = MAT2CELL(F) assumes is C a vector with all entries equal to 1 whose
%   length is equal to the number of components of F.
%
%   G = MAT2CELL(F, M, N) is similar to above, but allows three input arguments
%   so as to be consistent with the built-in MAT2CELL function.  If F is a
%   column CHEBFUN, then N takes the role of C above, and M should be the
%   scalar value 1.  If F is a row CHEBFUN, then M takes the role of C above,
%   and N should be the scalar value 1.
%
% Example:
%   f = chebfun(@(x) [sin(x), cos(x), exp(x), x], [0, pi])
%   g = mat2cell(f, 1, [1, 2, 1])
%
% See also NUM2CELL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return an empty result for empty inputs:
if ( isempty(f) )
    g = {f};
    return
end

% Get the size of the CHEBFUN array:
numFuns = numel(f.funs);
numComponents = size(f.funs{1}, 2);

% Parse the inputs:
if ( nargin == 1 )
    if ( f.isTransposed )
        M = ones(1, numComponents);
        N = 1;
    else
        M = 1;
        N = ones(1, numComponents);
    end
elseif ( nargin == 2 )
    if ( f.isTransposed )
        N = 1;
    else
        N = M;
        M = 1;
    end
end

% Handle row CHEBFUNs by transposing first.
if ( f.isTransposed )
    g = mat2cell(f.', N, M).'; % Note that M and N inputs get reversed.

    for k = 1:numel(g)
        g{k} = g{k}.';
    end

    return
end

% From this point on, f will be a column CHEBFUN.

% Check dimensions:
if ( ~isscalar(M) || (M ~= 1) || (sum(N) ~= numComponents) )
    error('CHEBFUN:CHEBFUN:mat2cell:size', ...
        ['Input arguments, M and N, must sum to each dimension of the', ...
        ' input size, [1,%d].'], numComponents);
end

% Call MAT2CELL of the FUNs:
cellFuns = cell(numFuns, numel(N));
for k = 1:numFuns
    cellFuns(k,:) = mat2cell(f.funs{k}, M, N);
end

% Create a cell which tells us which components are grouped together:
index = mat2cell(1:size(f.funs{1}, 2), M, N);

% Create a new CHEBFUN from each column (or row) of FUNs;
g = cell(1, numel(N));
for k = 1:numel(N)
    % Make the CHEBFUN.
    g{k} = chebfun(cellFuns(:,k));

    % Copy over higher-order impulses.
    g{k}.impulses = f.impulses(:,index{k},:);
    g{k} = tidyImpulses(g{k});
end

end
