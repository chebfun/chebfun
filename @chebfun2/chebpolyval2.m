function varargout = chebpolyval2( f, varargin )
%CHEBPOLYVAL2      Values on a tensor product grid.
%   X = CHEBPOLYVAL2(F) returns the matrix of values of F on a tensor
%   product grid.
%
%   [U, D, V] = CHEBPOLYVAL2(F) returns the low rank representation of the
%   values of F on a tensor product grid. X = U * D * V'.
%
%   [U, D, V] = CHEBPOLYVAL2(F,M,N) returns the values of F on a M-by-N
%   tensor product grid.
%
% See also CHEBCOEFFS2, PLOTCOEFFS2. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check. 
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

if ( nargin == 1 ) 
    % Get degrees:
    [m, n] = length( f );  
elseif ( nargin == 2 ) 
    error('CHEBFUN:CHEBFUN2:chebpolyval2:inputs', 'Dimension not specified.'); 
else
    m = varargin{ 1 }; 
    n = varargin{ 2 }; 
end

% Get the low rank representation for f. 
[cols, d, rows] = cdr(f);

tech = chebfunpref().tech(); 

C = tech.coeffs2vals(chebcoeffs( cols, n )); 
R = tech.coeffs2vals(chebcoeffs( rows, m )); 

% Evaluate: 
if ( nargout <= 1 )
    varargout = {C * d * R.'}; 
else
    varargout = {C , d, R}; 
end
    
end
