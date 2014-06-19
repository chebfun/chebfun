function r = rank(N)
%RANK   Rank of the partial differential operator.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

A = N.coeffs; 
if ( iscell(A) )
    error('CHEBFUN:CHEBOP2:rank:variableCoeffs', ...
        ['The operator has non-scalar variable coefficients. We do not ' ...
         'support this.'])
else
    r = rank(A);
end

end
