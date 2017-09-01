function varargout = qr(f, ignored)
%QR Orthogonal-triangular decomposition of a SEPARABLEAPPROX. 
% 
% [Q, R] = QR( F ), where F is a separableApprox, produces an unitary column
% quasimatrix Q and a upper-triangular row quasimatrix R so that F = Q * R. This
% is computed by a continuous analogue of QR. 
%
% [Q, R] = QR( F, 0 ) is the same as QR( F ). 
%
% [Q, R, E] = QR( F ) and [Q, R, E] = QR( F, 'vector') produces a vector E that 
% stores the pivoting locations. 
%
% For more information about this decomposition: 
% A. Townsend and L. N. Trefethen, Continuous analogues of matrix
% factorizations, Proc. Royal Soc. A., 2015. 
%
% See also LU, and CHOL. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) )
    varargout = cell(1, nargout); 
    return
end

% As always start with the CDR decomposition: 
[C, D, R] = cdr( f ); 

% Balance out the scaling, becareful about signs: 
sgns = sign( diag( D ) ); 
C = C * diag(sgns) * sqrt( abs(D) );   % Put the signs into Q
R = R * sqrt( abs(D) ); 

% QR of the column 
[Q, RC] = qr( C ); 

% Form R:  
R = RC * R.';   % still an upper-triangular quasimatrix! :)

% Output to user: 
if ( nargout <= 1 ) 
    varargout = { Q }; 
elseif ( nargout == 2 ) 
    varargout = { Q, R };
elseif ( nargout == 3 )
    E = f.pivotLocations;   % complete pivot Locations from GE. 
    varargout = { Q, R, E( :, 1 ) };  % Only pass the ones in x. 
end

end
