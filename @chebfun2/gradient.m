function F = gradient( f ) 
% GRADIENT   gradient of a chebfun2

fx = diff(f, 1, 2);   % diff in x-variable
fy = diff(f, 1, 1);   % diff in y-variable 

F = chebfun2v( {fx, fy} );

end 