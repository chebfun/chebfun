function F = getVals2CoeffsTransform(~)
%GETVALS2COEFFSTRANSFORM   Returns the values to coeffs transform on the sphere.
%
% See also SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = @(u) reshape(trigtech.vals2coeffs(trigtech.vals2coeffs(u).').', ...
    length(u)^2, 1);

end