function r = rank(f, tol)
%RANK  Rank of a CHEBFUN2.
%   RANK(F) produces an estimate of the rank of the approximant F.
%
%   RANK(F, TOL) is the number of singular values of F greater than TOL/N, where
%   N is the first singular value of F.
%
% See also LENGTH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( f ) )
    r = [];
    return
end

% HARD tolerance. 
if ( nargin == 1 )
    tol = 0;
end

% Compute the singular values of f. 
s = svd( f ); 

% Check for zero function.
if ( max( s ) == 0  )  
    r = 0; 
else
    % r = no. of s.v. above relative tol.
    r = find( s/s(1) > tol, 1, 'last');  
end

end
