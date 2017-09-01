function D = diff(disc, m)
%DIFF   Differentiation operator for TRIGSPEC.
%   D = DIFF(DISC, m) returns a differentiation operator for the TRIGSPEC
%   class that represents the mth derivative.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Useful information:
dom = disc.domain;
N = disc.dimension;

if ( m == 0 )          
    % 0th order derivative is easy!
    D = speye(sum(N));
    
else
    % Rescale the differentiation matrix.
    rescaleFactor = dom(2) - dom(1);
    D = trigspec.diffmat(N, m) * (2*pi/rescaleFactor)^m;
end

end
