function varargout = lu( f, thresh )
%LU   LU factorization of a SEPARABLEAPPROX.
%
% [L, U] = LU( F ) returns two quasimatrices L and U of size inf by k and 
% k by inf, respectively, where k is the rank of the separableApprox F. 
% The quasimatrices L and U are "psychologically" lower and upper triangular.
% L is also unit lower triangular. This is computed by a continuous analogue of
% Gaussian elimination with complete pivoting. 
%
% [L, U, P] = lu( F ) returns a k by 2 matrix P, containing the pivoting elements.
%
% [L, U, P, Q] = LU( F ) returns two vectors P and Q containing the x-values and
% y-values of the pivoting elements.
%
% [L, U, P] = LU(F, THRESH) returns the LU factorization of F that removes 
% the tail of pivots below THRESH.
% 
% For more information about the factorization: 
% A. Townsend and L. N. Trefethen, Continuous analogues of matrix
% factorizations, submitted, 2014. 
%
% See also CHOL, QR. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    thresh = 0; 
end

% Empty check: 
if ( isempty( f ) ) 
    varargout = cell(1, nargin);   % empty output 
    return
end

% Get the CDR factorization
[C, D, R] = cdr( f ); 

% Pivot locations 
PivLoc = f.pivotLocations;

% Pivot sizes 
PivSize = f.pivotValues;

% rank of f. 
k = length( f ); 

% Make C unit lower triangular 
Scl = diag( C( PivLoc(:, 2), : ) );
C = C * spdiags( 1./Scl, 0, k, k );
R = R * ( spdiags( Scl, 0, k, k ) * D );

L = C; 
U = R; 

% chop the LU: 
if ( thresh > 0 )
   newk = find( abs( PivSize ) < thresh, 1, 'first' );
   L = L( :, 1 : newk );
   U = U( :, 1 : newk );
end

% Convert to row quasimatrix: 
U = U.'; 

% Output to the user: 
if ( nargout < 2 ) 
    varargout = { L }; 
elseif ( nargout == 2 )
    varargout = { L, U }; 
elseif ( nargout == 3 ) 
    varargout = { L, U, PivLoc };
elseif ( nargout == 4 )
    varargout = { L, U, PivLoc(:,1), PivLoc(:,2) };
end

end
