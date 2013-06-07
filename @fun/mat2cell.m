function g = mat2cell(f, M, N)
%MAT2CELL  Convert an array-valued FUN into an ARRAY of FUN objects.
%
%   G = MAT2CELL(F, C) breaks up the array-valued FUN F into a single row
%   cell array G of FUN objects. C is the vector of column sizes, and if F
%   has COL columns this must sum to COL. The elements of C determine the size
%   of each cell in G, subject to the following formula for I = 1:LENGTH(C),
%               SIZE(C{I},2) == C(I).
%
%   G = MAT2CELL(F) will assume that C = ones(1, COL).
%
%   G = MAT2CELL(F, M, N) is similar to above, but allows three input arguments
%   so as to be consistent with the built in MAT2CELL function. Here N takes the
%   role of C above, and M should be the scalar value 1.
%
% Example:
%   f = bndfun(@(x) [sin(x), cos(x), exp(x), x], [0 pi])
%   g = mat2cell(f, 1, [1 2 1])
%
% See also CELL2MAT.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return an empty result:
if ( isempty(f) )
    g = f;
    return
end

% Get the size of the values matrix:
[ignored, m] = size(f.onefun);

% Parse the inputs:
if ( nargin == 1 )
    M = 1;
    N = ones(1, m);
elseif ( nargin == 2 ) 
    N = M;
    M = 1;
end

if ( ~isscalar(M) || M ~= 1)
    error('CHEBFUN:FUN:mat2cell:size', ...
        ['Input arguments, D1 through D2, must sum to each dimension of the', ...
        ' input matrix size, [1,%d].'], size(f.values, 2));
end

% Split the values and the coeffs into cells of the correct size:
oneFuns = mat2cell(f.onefun,M,N);

% Create a cell for storing the FUN objects
g = cell(1, numel(N));

% Append the data to the new entries in the cell:
for k = 1:numel(N)
    % Create a FUN
    gk = f.make();
    
    % Assign values to the fields of the FUN
    gk.domain = f.domain;
    gk.onefun = oneFuns{k};
    gk.mapping = f.mapping;

    % Store the CHEBTECH in the cell-array returned.
    g{k} = gk;
end

end