function out = vscale(f)
%VSCALE    Vertical scale of a CHEBFUN object.
%   VSCALE(F) returns an estimate of the maximum absolute value of F. VSCALE
%   always returns a scalar, even when F is an array-valued CHEBFUN. Vertical
%   scales of each of the piecewise components and columns of F are given by
%   get(F, 'vscale-local');
%
% See also MAX.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the local vscales:
v = get(f, 'vscale-local');

% Compute the maximum:
out = max(v(:));

end