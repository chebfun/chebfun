function varargout = svd( f )
%SVD of a chebfun2.  
%
% SVD(F) returns the singular values of F. The number of singular values
% returned is equal to the rank of the chebfun2. 
%
% S = SVD(F) returns the singular values of F. S is a vector of singular
% values in decreasing order. 
% 
% [U S V] = SVD(F) returns the SVD of F. U and V are quasi-matrices of 
% orthogonal chebfuns and S is a diagonal matrix with the singular values
% on the diagonal.
%
% The length and rank of a chebfun2 are slightly different quantities. 
% LENGTH(F) is the number of pivots used by the Chebfun2 constructor, and 
% RANK(F) is the number of significant singular values of F. The relation 
% RANK(F) <= LENGTH(F) should always hold. 


if isempty(F) % check for empty chebfun2. 
    varargout = {chebfun,[],chebfun};
    return
end

dom = f.domain; 
width = diff( dom( 1:2 ) ); 
height = diff( dom( 3:4 ) ); 
pivots = f.pivotValues;

% If the function is the zero function then special care is required. 
if ( norm( pivots ) == 0 )
    if ( nargout > 1 ) 
        U = chebfun2( 1./sqrt( width ), dom(1:2) ); 
        V = chebfun2( 1./sqrt( height ), dom(3:4) );
        varargout = { U, 0, V }; 
    else
        varargout = { 0 };
    end
else

% If the function is non-zero then do the standard stuff. 
[Qleft, Rleft] = qr( f.cols ); 
[Qright, Rright] = qr( f.rows );
[U, S, V] = svd( Rleft * diag( 1./f.pivotValues ) * Rright.' );
U = Qleft * U; 
V = Qright * V; 

% Output just like the svd of a matrix. 
if ( nargout > 1 )
    varargout = { U, S, V};
else
    varargout = { S };
end

end