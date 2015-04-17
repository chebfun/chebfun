function f = uminus( f ) 
% UMINUS      Unary minus for a spherefun 
% 
%  F = UMINUS( F )   

f.BlockDiag = -f.BlockDiag; 

end