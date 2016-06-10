function out = ndf(f)
%NDF   Number of degrees of freedom (parameters) needed to represent a 
%   CHEBFUN3T.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    out = 0;
else
    out = numel(f.coeffs);
end

end