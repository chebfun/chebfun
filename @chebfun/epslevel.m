function out = epslevel(F, ignoreUnhappy) %#ok<INUSD>
%EPSLEVEL   Accuracy estimate of a CHEBFUN object.
%   EPSLEVEL(F) returns an estimate of the relative error in the CHEBFUN F. This
%   is defined as the maximum of the product of the local vscales and epslevels,
%   divided by the global vscale.
%
%   EPSLEVEL(F, 'ignoreUnhappy') ignores local epslevels of unhappy subintervals
%   when computing the global estimate.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = eps;

end
