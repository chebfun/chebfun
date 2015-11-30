function varargout = svd( f )
%SVD    Singular value decomposition of a SEPARABLEAPPROX.
%   SVD(F) returns the singular values of F. The number of singular values
%   returned is equal to the rank of the SEPARABLEAPPROX.
%
%   S = SVD(F) returns the singular values of F. S is a vector of singular
%   values in decreasing order.
%
%   [U, S, V] = SVD(F) returns the SVD of F. U and V are quasi-matrices of
%   orthogonal CHEBFUN objects and S is a diagonal matrix with the singular
%   values on the diagonal.
%
%   The length and rank of a SEPARABLEAPPROX are slightly different quantities.
%   LENGTH(F) is the number of pivots used by the constructor, and
%   RANK(F) is the number of significant singular values of F. The relation
%   RANK(F) <= LENGTH(F) should always hold.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) )
    varargout = { chebfun, [], chebfun };
    return
end

% if ( iszero( f ) ) 
%     varargout = {0}; 
%     return
% end

% Get the low rank representation for f.
[cols, D, rows] = cdr(f);
d = diag(D);

% Extract information:
dom = f.domain;
width = diff( dom( 1:2 ) );
height = diff( dom( 3:4 ) );

% If the function is the zero function then special care is required.
if ( norm( d ) == 0 )
    if ( nargout > 1 )
        f = 1 + 0*f;
        U = 1/sqrt( width )*simplify(f.cols, [], 'globaltol');
        V = 1/sqrt( height )*simplify(f.rows, [], 'globaltol');
        varargout = { U, 0, V };
    else
        varargout = { 0 };
    end
    
else
    
    % If the function is non-zero then do the standard stuff.
    %
    % Algorithm:
    %   f = C D R'                 (cdr decomposition)
    %   C = Q_C R_C                (qr decomposition)
    %   R = Q_R R_R                (qr decomposition)
    %   f = Q_C (R_C D R_R') Q_R'
    %   R_C D R_R' = U S V'        (svd)
    
    [Qleft, Rleft] = qr( cols );
    [Qright, Rright] = qr( rows );
    [U, S, V] = svd( Rleft * D * Rright.' );
    U = Qleft * U;
    V = Qright * V;
    
    % Output just like the svd of a matrix.
    if ( nargout > 1 )
        varargout = { U, S, V };
    else
        varargout = { diag( S ) };
    end
    
end

end
