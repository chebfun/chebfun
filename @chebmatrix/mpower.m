function B = mpower(A, pow)
%^   Repeated composition of a CHEBMATRIX.
%
% See also MTIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( (pow ~= round(pow)) || (pow < 0) )
    error('CHEBFUN:CHEBMATRIX:power:badPower', ...
        'Power must be a positive integer.')
end

% Create an "identity" CHEBMATRIX for the given variable types so we can start
% the repeated composition.
B = identity(A);
       
for i = 1:pow
    B = B * A;
end

end
