function r = rank(N)
%RANK   Rank of the partial differential operator.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

A = N.coeffs; 
if ( iscell(A) )
    if ( ~isempty(N.S) ) 
        r = size(N.S, 1); 
    else
        % attempt to compute it: 
        [U, S, V] = chebop2.separableFormat(N); 
        r = size(S, 1); 
    end
else
    r = rank(A);
end

end
