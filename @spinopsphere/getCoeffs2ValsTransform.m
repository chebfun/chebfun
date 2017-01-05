function F = getCoeffs2ValsTransform(~)
%GETCOEFFS2VALSTRANSFORM   Returns the coeffs to values transform on the sphere.
%
% See also SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = @(u) trigtech.coeffs2vals(trigtech.coeffs2vals( ...
    reshape(u, sqrt(length(u)), sqrt(length(u)))).').';

end