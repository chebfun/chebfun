function T = spconvert(n, lam)
%SPCONVERT   Compute sparse representation for conversion operators. 
%   CONVERMAT(N, LAM) returns the truncation of the operator that transforms
%   C^{lam} (Ultraspherical polynomials) to C^{lam+1}.  The truncation gives
%   back a matrix of size n x n.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Relation is: C_n^(lam) = (lam/(n+lam))(C_n^(lam+1) - C_{n-2}^(lam+1))

if ( lam == 0 )
    dg = .5*ones(n - 2, 1);
    T = spdiags([1 0 ; .5 0 ; dg -dg], [0 2], n, n);
else
    dg = lam./(lam + (2 : n - 1))';
    T = spdiags([1 0 ; lam./(lam + 1) 0 ; dg -dg], [0 2], n, n);
end

end
