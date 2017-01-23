function varargout = chebpolyval3(f, varargin)
%CHEBPOLYVAL3   Values of a CHEBFUN3 object F on a tensor product grid.
%   X = CHEBPOLYVAL3(F, M, N, P) returns an M x N x P tensor of values of a
%   CHEBFUN3 object F on a tensor product grid.
%
%   X = CHEBPOLYVAL3(F) does the same but sets M, N, and P equal to the 
%   length of F.
%
%   [CORE, C, R, T] = CHEBPOLYVAL3(F) returns the low rank representation 
%   of the values of F on a tensor product grid, i.e.,
%   X = CORE x_1 C x_2 R x_3 T.
% 
%   Example: For an order-3 discrete tensor R, we should have 
%                                       R \approx chebpolyval3(chebfun3(R)).

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check. 
if ( isempty(f) )
    varargout = { [] }; 
    return
end

if ( nargin == 1 ) 
    % Get degrees:
    [m, n, p] = length(f);
elseif ( nargin ~= 4 ) 
    error('CHEBFUN:CHEBFUN3:chebpolyval3:inputs', 'Dimension not specified.'); 
else
    m = varargin{1}; 
    n = varargin{2}; 
    p = varargin{3}; 
end

% Get low rank representation of f:
[core, cols, rows, tubes] = tucker(f);

tech = chebfunpref().tech(); 
colVals = tech.coeffs2vals(chebcoeffs(cols, m));
rowVals = tech.coeffs2vals(chebcoeffs(rows, n));
tubeVals = tech.coeffs2vals(chebcoeffs(tubes, p));

% Evaluate: 
if ( nargout <= 1 )
    varargout = {chebfun3.txm(chebfun3.txm(chebfun3.txm(core, colVals, 1), ...
        rowVals, 2), tubeVals, 3)};
else
    varargout = {core, colVals, rowVals, tubeVals};
end

end