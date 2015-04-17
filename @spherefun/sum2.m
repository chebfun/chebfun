function v = sum2( f ) 
% Definite integration of a spherefun. 

cs = cumsum( f.Cols );
v = feval(cs, 0)  * f.BlockDiag * sum(f.Rows).';

end