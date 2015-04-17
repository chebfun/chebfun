function [C, D, R] = cdr( f ) 
% CDR     CDR decomposition of a spherefun 
% 
%   [C, D, R] = CDR( F ) returns the CDR decomposition of a spherefun. 

C = f.Cols; 
R = f.Rows; 
D = f.BlockDiag; 

end