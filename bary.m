function fx = bary(x, fvals, xk, vk)
%BARY   Barycentric interpolation formula.
%   BARY(X, FVALS, XK, VK) uses the 2nd form barycentric formula with weights VK
%   to evaluate an interpolant of the data {XK, FVALS(:,k)} at the points X.
%   Note that X, XK, and VK should be column vectors, and FVALS, XK, and VK
%   should have the same length.
%
%   BARY(X, FVALS) assumes XK are the 2nd-kind Chebyshev points and VK are the
%   corresponding barycentric weights.
%
% See also CHEBTECH.CLENSHAW.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: This is simply a wrapper for CHEBTECH.BARY().
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs:
n = size(fvals, 1);

% Default to Chebyshev nodes and barycentric weights:
if ( nargin < 3 )
    xk = chebtech2.chebpts(n);
end
if ( nargin < 4 )
    vk = chebtech2.barywts(n);
end
    
% Call CHEBETCH.BARY.
fx = chebtech.bary(x, fvals, xk, vk);

end
