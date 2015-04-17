function varargout = svd( f ) 
% SVD      Singular value decomposition of a spherefun 
% 
% S = SVD( F )  returns the singular values of F 
% 
% [U, S, V] = svd( F ) returns the singular value decomposition of F. 

% COMPRESS is actually doing the SVD: 
f = compress( f ) ; 

if ( nargout <= 1 ) 
    varargout = { full( diag( f.BlockDiag ) ) }; 
elseif ( nargout == 3 ) 
    varargout = { f.Cols, f.BlockDiag, f.Rows }; 
end

end