function normF = norm(f, p)
%NORM       Frobenius norm of a CHEBFUN3T object.
%    NORM(F) = NORM(F,'fro') = sqrt(triple integral of abs(F)^2).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    % Default to the Frobenius norm.
    p = 'fro';
end

if ( isempty(f) )
    % Empty chebfun3T has no norm.
    normF = [];
elseif strcmp(p, 'fro')
    normF = sqrt(sum3(f.^2));  
else
    error('CHEBFUN:CHEBFUN3T:norm:norm', ...
        'CHEBFUN3T does not support this norm.');
end

end