function G = mat2cell(F, M, N)
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

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return an empty result for empty inputs:
if ( isempty(F) )
    G = {F};
    return
end

% Get the size of the CHEBFUN array:
numCols = numColumns(F);

% Parse the inputs:
if ( nargin == 1 )
    if ( F(1).isTransposed )
        M = ones(1, numCols);
        N = 1;
    else
        M = 1;
        N = ones(1, numCols);
    end
elseif ( nargin == 2 )
    if ( F(1).isTransposed )
        N = 1;
    else
        N = M;
        M = 1;
    end
end

% Handle row CHEBFUNs by transposing first.
if ( F(1).isTransposed )
    G = mat2cell(F.', N, M).'; % Note that M and N inputs get reversed.

    for k = 1:numel(G)
        G{k} = G{k}.';
    end

    return
end

% From this point on, f will be a column CHEBFUN.

% Check dimensions:
if ( ~isscalar(M) || (M ~= 1) || (sum(N) ~= numCols) )
    error('CHEBFUN:CHEBFUN:mat2cell:size', ...
        ['Input arguments, M and N, must sum to each dimension of the', ...
        ' input size, [1,%d].'], numCols);
end

if ( numel(F) > 1 )
    G = quasiMat2cell(F, M, N);
else
    G = columnMat2cell(F, M, N);
end
    
end

function g = columnMat2cell(f, M, N)
% MAT2CELL for an array-valued CHEBFUN.

% Call MAT2CELL of the FUNs:
numFuns = numel(f.funs);
cellFuns = cell(numFuns, numel(N));
for k = 1:numFuns
    cellFuns(k,:) = mat2cell(f.funs{k}, M, N);
end

% Create a cell which tells us which components are grouped together:
index = mat2cell(1:size(f.funs{1}, 2), M, N);

% Create a new CHEBFUN from each column (or row) of FUNs and store it in a CELL;
for k = numel(N):-1:1
    % Make the CHEBFUN.
    g{k} = chebfun(cellFuns(:,k).');

    % Copy over pointValues.
    g{k}.pointValues = f.pointValues(:,index{k});    
end

end

function g = quasiMat2cell(f, M, N)
% MAT2CELL for a CHEBFUN array.

% Get the correct indices fron the built-in MAT2CELL():
idx = mat2cell(1:numel(f), M, N);
% Store these coluns in a CELL array:
for k = numel(idx):-1:1
    g{k} = f(idx{k});
end

end
