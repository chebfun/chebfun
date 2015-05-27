function v = sum2( f ) 
% Definite integration of a spherefun. 

meanCols = sum( f.cols )/2;  % mean(f.cols) 
cs = pi*cumsum( f.cols );
v = (feval(cs, 0)+2*pi*meanCols) * f.blockDiag * sum(f.rows).';

end