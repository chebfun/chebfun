function vscl = vscale( f )
%VSCALE   Estimate the vertical scale of a function. 
%   VSCALE(F) estimates the vertical scale (also known as the dynamical range)
%   of a function. This is required because a CHEBTECH does not store its
%   interpolation data. If F is an array-valued CHEBTECH with K columns, then
%   the result is a row vector of length K.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the coefficients:
c = f.coeffs;

% Check isempty: 
if ( isempty( c ) ) 
    vscl = 0;
elseif ( size(c, 1) == 1 )
    vscl = abs(c);
else
    % Compute values:
    vals = f.coeffs2vals( c );
    % Take max:
    vscl = max(abs(vals), [], 1);
end

end
