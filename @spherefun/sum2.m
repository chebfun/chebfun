function v = sum2( f ) 
% Definite integration of a spherefun. 

meanCols = sum( f.Cols )/2;  % mean(f.Cols) 
cs = pi*cumsum( f.Cols );
v = (feval(cs, 0)+2*pi*meanCols) * f.BlockDiag * sum(f.Rows).';

end