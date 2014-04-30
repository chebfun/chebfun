function out = vscale(F)
%VSCALE   Vertical scale of a CHEBFUN object.
%   VSCALE(F) returns an estimate of the maximum absolute value of F. VSCALE
%   always returns a scalar, even when F is an array-valued CHEBFUN or a
%   quasimatrix. Vertical scales of each of the piecewise components and columns
%   of F are given by get(F, 'vscale-local');
%
% See also MAX, MINANADMAX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) )
    out = [];
    return
end

out = 0;
for k = 1:numel(F)
    % Get the local vscales:
    v = get(F(k), 'vscale-local');

    % Compute the maximum:
    out = max(out, max(v(:)));
end

end
