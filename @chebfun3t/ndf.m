function out = ndf(f)
%NDF   Number of degrees of freedom needed to represent a CHEBFUN3T object.
%
% See also CHEBFUN3/NDF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    out = 0;
else
    out = numel(f.coeffs);
end

end