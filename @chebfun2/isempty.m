function bol = isempty( f ) 
% ISEMPTY True for empty chebfun2

bol_cols = isempty( f.cols ); 
bol_rows = isempty( f.rows ); 

bol = bol_cols & bol_rows; 

end