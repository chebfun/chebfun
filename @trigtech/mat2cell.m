function g = mat2cell(f, M, N)
%MAT2CELL   Convert an array-valued TRIGTECH into a cell array of TRIGTECH
%objects.
%
%   G = MAT2CELL(F, C) breaks up the array-valued TRIGTECH F into a single row
%   cell array G of TRIGTECH objects. C is the vector of column sizes, and if F
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
%   f = trigtech(@(x) [sin(x), cos(x), exp(cos(x)), sin(sin(x)])
%   g = mat2cell(f, 1, [1 2 1])
%
% See also CELL2MAT.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Return an empty result:
if ( isempty(f) )
    g = f;
    return
end

% Get the size of the values matrix:
[n, m] = size(f.values);

% Parse the inputs:
if ( nargin == 1 )
    M = 1;
    N = ones(1, m);
elseif ( nargin == 2) 
    N = M;
    M = 1;
end

if ( ~isscalar(M) || (M ~= 1) )
    error('CHEBFUN:TRIGTECH:mat2cell:size', ...
        ['Input arguments M and N must sum to each dimension of the', ...
        ' input matrix size, [1,%d].'], size(f.values, 2));
end

% Split the values and the coefficients into cells of the correct size:
values = mat2cell(f.values, n, N);
coeffs = mat2cell(f.coeffs, n, N);
vscale = mat2cell(f.vscale, 1, N);
epslevel = mat2cell(f.epslevel, 1, N);
isReal = mat2cell(f.isReal, 1, N);

% Create a cell for storing the TRIGTECH objects
g = cell(1, numel(N));

% Append the data to the new entries in the cell:
for k = 1:numel(N)
    % Create a TRIGTECH:
    gk = f.make();
    
    % Assign values to the fields of the TRIGTECH:
    gk.ishappy = f.ishappy;
    gk.values = values{k};
    gk.coeffs = coeffs{k};
    gk.vscale = vscale{k};
    gk.epslevel = epslevel{k};
    gk.isReal = isReal{k};
    
    % Store the TRIGTECH in the cell-array returned.
    g{k} = gk;
end

end