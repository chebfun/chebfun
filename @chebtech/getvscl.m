function vscl = getvscl( f )
%GETVSCL(F)   Estimate the vertical scale of a function. 
%   VSCALE = GETVSCL(F) estimates the vertical scale (also known as the
%   dynamical range) of a function. This is required because a CHEBTECH does not
%   store its interpolation data. If F is an array-valued CHEBTECH with K
%   columns, then VSCALE is a row vector of length K.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check isempty: 
if ( isempty( f.coeffs ) ) 
    vscl = 0;
else
    % compute values
    values = f.coeffs2vals( f.coeffs );
    vscl = max(abs(values), [], 1);
end

end