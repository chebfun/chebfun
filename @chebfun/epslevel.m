function out = epslevel(F, ignoreUnhappy) %#ok<INUSD>
%EPSLEVEL   Accuracy estimate of a CHEBFUN object.
%   EPSLEVEL(F) returns an estimate of the relative error in the CHEBFUN F. This
%   is defined as the maximum of the product of the local vscales and epslevels,
%   divided by the global vscale.
%
%   EPSLEVEL(F, 'ignoreUnhappy') ignores local epslevels of unhappy subintervals
%   when computing the global estimate.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    ignoreUnhappy = 0;
else
    ignoreUnhappy = 1;
end

out = 0;
for k = 1:numel(F)
    out = max(out, columnEpslevel(F(k), ignoreUnhappy));
end

% Make it a _relative_ error estimate:
vs = vscale(F);
vs(vs < out) = 1;
out = out / vs;

end

function out = columnEpslevel(f, ignoreUnhappy)

% Get the local vscales:
v = get(f, 'vscale-local');
% Get the local epslevels:
e = get(f, 'epslevel-local');
% Ignore unhappy intervals?
if ( ignoreUnhappy )
    ish = cellfun(@(f) get(f, 'ishappy'), f.funs);
    e(~ish) = 0;
end
% Compute their product:
ve = bsxfun(@times, v, e);
% Compute the maximum of their product:
out = max(ve(:));

% [TODO]: Remove this hack!
if ( isnan(out) || ~logical(out) )
    out = chebfunpref().eps;
end

end
