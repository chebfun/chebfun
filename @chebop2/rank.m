function r = rank(N)
% RANK  rank of the partial differential operator.

A = N.coeffs; 
if iscell(A)
   error('The operator has non-scalar variable coefficients. We do not support this.')
else
    r = rank( A );
end
end