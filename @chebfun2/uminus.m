function F = uminus( F )
%UMINUS   Unary minus for a CHEBFUN2. 

F.pivotValues = -F.pivotValues;

end