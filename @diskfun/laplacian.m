function f = laplacian( f ) 
% LAPLACIAN     Scalar laplacian of a diskfun 

f = diff(f, 1, 2) + diff( f, 2, 2); 

end 