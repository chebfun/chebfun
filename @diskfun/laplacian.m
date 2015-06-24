function f = laplacian( f ) 
% LAPLACIAN     Scalar laplacian of a spherefun 
% 
% F = LAPLACIAN( F ) 

f = diff(f, 1, 2) + diff( f, 1, 1); 

end 