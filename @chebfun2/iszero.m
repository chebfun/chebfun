function bol = iszero( f ) 
% ISZERO   Check if a chebfun2 is identically zero on its domain. 

pivots = f.pivotValues; 
cols = f.cols; 
rows = f.rows; 

% Check pivot values
bol = ( norm( pivots, inf ) == 0 );

rk = length( f );
if ( bol ) 
    bolcols = zeros( length( f ), 1 );
    bolrows = zeros( length( f ), 1 );
    for j = 1 : rk
        bolcols( j ) = iszero( cols(:, j) ); 
        bolrows( j ) = iszero( rows(:, j) ); 
    end
    
    bolslices = ( any(bolcols) & any(bolrows) );
    
    bol = bolslices & bol; 
end



end