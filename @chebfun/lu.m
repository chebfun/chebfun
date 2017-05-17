function [L, U, p] = lu( A , varargin )
%LU  LU factorization of a quasimatrix. 
% 
% [L, U] = lu(A) stores an upper-triangular matrix U and a "psychologically"
% lower triangular quasimatrix in L so that A = L*U. Here, A is a column 
% quasimatrix. 
%
% [L, U, p] = lu(A) returns the pivoting information so that L(p,:) is a lower
% triangular matrix. 
%
% For more information about the factorization: 
% A. Townsend and L. N. Trefethen, Continuous analogues of matrix
% factorizations, submitted, 2014. 
%
% See also QR. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( A ) ) 
    L = chebfun; 
    U = []; 
    p = []; 
    return
end

if ( A(1).isTransposed )
    error('CHEBFUN:CHEBFUN:lu:sizes',...
        'CHEBFUN LU works only for column CHEBFUN objects.');
end

L = A; 
U = zeros( size(L, 2) ); 
p = NaN( size(L, 2), 1);
% Do GE with partial pivoting: 
for j = 1 : size(A, 2)
    Acol = extractColumns(A, j);
    
    % find maximum absolute entry: 
    [vals, pos] = minandmax( Acol );
    [ignored, idx] = max( abs( vals ) ); 
    pos = pos( idx );
    mx = feval( Acol, pos ); 
    Arow = feval(A, pos);     
    
    if ( ismember(pos, p) )
        error('CHEBFUN:CHEBFUN:lu:pivot',...
            'Duplicated pivot location, likely due to ill-conditioning.');
    end
    
    % Store upper-triangular part: 
    U( j, : ) = Arow; 
    if ( j == 1 ) 
        L = Acol / mx; 
    else
        L = horzcat(L, Acol / mx ); 
    end
    
    % Do one step of GE: 
    A = A - Acol * Arow / mx; 
    
    % store pivot locations: 
    p( j ) = pos;
end

% Make upper-triangular even with rounding errors:
U = triu( U );

end
