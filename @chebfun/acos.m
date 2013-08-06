function g = acos(f, pref)
%ACOS   Inverse cosine of a chebfun.
%
% See also COS, ACOSD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @acos, pref);

end