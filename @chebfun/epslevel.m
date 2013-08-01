function out = epslevel(f)
%EPSLEVEL    Accuracy estimate of a CHEBFUN object.
%   EPSLEVEL(F) returns an estimate of the relative error CHEBFUN F. This is
%   defined as the maximum of the product of the local vscales and epslevels,
%   i.e.,  max(max(get(F, 'vscale-local').*get(f, 'epslevel-local')));

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the local vscales:
v = get(f, 'vscale-local');
% Get the local epslevels:
e = get(f, 'epslevel-local');
% Compute their product:
ve = bsxfun(@times, v, e);
% Compute the maximum of their product:
out = max(ve(:));
% Make is relative error:
out = out / vscale(f);

end