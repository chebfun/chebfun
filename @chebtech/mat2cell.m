function g = mat2cell(f, M, N)
%MAT2CELL   Convert an array-valued CHEBTECH into a cell array of CHEBTECH
%           objects.
%   G = MAT2CELL(F, C) breaks up the array-valued CHEBTECH F into a single row
%   cell array G of CHEBTECH objects. C is the vector of column sizes, and if F
%   has COL columns this must sum to COL. The elements of C determine the size
%   of each cell in G, subject to the following formula for I = 1:LENGTH(C),
%   SIZE(G{I},2) == C(I).
%
%   G = MAT2CELL(F) will assume that C = ones(1, COL).
%
%   G = MAT2CELL(F, M, N) is similar to above, but allows three input arguments
%   so as to be consistent with the built-in MAT2CELL function. Here N takes the
%   role of C above, and M should be the scalar value 1.
%
% Example:
%   f = chebtech.constructor(@(x) [sin(x), cos(x), exp(x), x])
%   g = mat2cell(f, 1, [1 2 1])
%
% See also CELL2MAT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return an empty result:
if ( isempty(f) )
    g = f;
    return
end

% Get the size of the values matrix:
[n, m] = size(f.coeffs);

% Parse the inputs:
if ( nargin == 1 )
    M = 1;
    N = ones(1, m);
elseif ( nargin == 2) 
    N = M;
    M = 1;
end

if ( ~isscalar(M) || (M ~= 1) )
    error('CHEBFUN:CHEBTECH:mat2cell:size', ...
        ['Input arguments M and N must sum to each dimension of the', ...
        ' input matrix size, [1,%d].'], size(f.coeffs, 2));
end

% Split the coefficients into cells of the correct size:
coeffs = mat2cell(f.coeffs, n, N);

% Create a cell for storing the CHEBTECH objects
g = cell(1, numel(N));

% Append the data to the new entries in the cell:
for k = 1:numel(N)
    % Create a CHEBTECH
    gk = f.make();
    
    % Assign values to the fields of the CHEBTECH
    gk.ishappy = f.ishappy;
    gk.coeffs = coeffs{k};
    
    % Store the CHEBTECH in the cell-array returned.
    g{k} = gk;
end


end
