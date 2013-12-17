function f = simplify( f ) 
% Simplify a chebfun2 

% Simplify the column and row slices. 
f.cols = simplifty( f.cols ); 
f.rows = simplifty( f.rows ); 

% Should this also simplify the rank?  More difficult. 

end 