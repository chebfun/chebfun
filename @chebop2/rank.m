function r = rank(N)
% RANK  rank of the PDE operator.

A = N.coeffs; 
if iscell(A)
   error('The operator has non-scalar variable coefficients. At the moment we do not support this.')
else
    r = rank(A);
end
end