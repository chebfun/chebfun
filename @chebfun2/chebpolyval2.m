function varargout = chebpolyval2( f, varargin )
%CHEBPOLYVAL2   Values on a tensor Chebyshev grid.
%   X = CHEBPOLYVAL2(F) returns the matrix of values of F on a Chebyshev tensor
%   grid.
%
%   [U, D, V] = CHEBPOLYVAL2(F) returns the low rank representation of the
%   values of F on a tensor Chebyshev grid. X = U * D * V'.
%
%   [U, D, V] = CHEBPOLYVAL2(F,M,N) returns the values of F on a M-by-N
%   Chebyshev tensor grid.
%
% See also CHEBPOLY2, CHEBPOLYPLOT2. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check. 
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

if ( nargin == 1 ) 
    % Get degrees:
    [m, n] = length( f );  
elseif ( nargin == 2 ) 
    error('CHEBFUN2:CHEBPOLYVAL2:INPUTS', 'Dimension not specified.'); 
else
    m = varargin{ 1 }; 
    n = varargin{ 2 }; 
end

% Get the low rank representation for f. 
[cols, d, rows] = cdr(f);

% Evaluate the columns / rows
C = resize( chebpoly( cols ).', n );
R = resize( chebpoly( rows ).', m );

% Convert these values to coefficients.
C = chebtech2.coeffs2vals( C );
R = chebtech2.coeffs2vals( R );

% Evaluate: 
if ( nargout <= 1 )
    varargout = {C * d * R.'}; 
else
    varargout = {C , d, R}; 
end
    

end

function X = resize( X, N )
% Resize the matrix to have length N.
    
% Get size:
[mX, nX] = size( X ); 

if ( mX > N ) 
    % Truncate:
    X = X((end-N+1):end, :); 
else
    % Pad:
    X = [zeros(N - mX, nX) ; X]; 
end

end