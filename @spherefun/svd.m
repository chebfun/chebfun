function varargout = svd( f ) 
% SVD      Singular value decomposition of a spherefun 
% 
% S = SVD( F )  returns the singular values of F 
% 
% [U, S, V] = svd( F ) returns the singular value decomposition of F. 

% TODO: Check if this is good enough, structure is probably not preserved.
% :( 

f = compress(f);

C = f.Cols; 
R = f.Rows; 
D = f.BlockDiag; 

% QR the outside quasimatrices: 
[QC, RC] = qr( C, 0 ); 
[QR, RR] = qr( R, 0 ); 

% Now do inner matrix: 
[U, S, V] = svd( RC * D * RR.' ); 

% Put it together: 
QC = QC * U; 
RC = QR * V;

if ( nargout <= 1 ) 
    varargout = { full( diag( S ) ) }; 
elseif ( nargout == 3 ) 
    varargout = { QC, S, RC }; 
end

end