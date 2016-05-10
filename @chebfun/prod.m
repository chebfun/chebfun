function out = prod(F, dim)
%PROD   Product integral.
%   PROD(F) is the definite product integral of the CHEBFUN F, which is defined
%   as exp(sum(log(F))).
%
%   PROD(F, DIM) takes products along the dimension DIM of the quasimatrix F.
%   If DIM is the continuous dimension of F, the result is the product integral
%   of each column or row of F, as described above.  Otherwise, it is the
%   pointwise product across the columns or rows of F.
%
%   PROD(F) and PROD(F, DIM) are exactly equivalent to EXP(SUM(LOG(F))) and
%   EXP(SUM(LOG(F), DIM)), respectively.
%
% See also CUMPROD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Determine the continuous dimension.
if ( F(1).isTransposed )
    continuousDim = 2;
else
    continuousDim = 1;
end

% Multiply in the continuous dimension by default.
if ( nargin < 2 )
    dim = continuousDim;
end

if ( dim == continuousDim )
    % Compute the product integrals.
    out = exp(sum(log(F), dim));
else
    % Multiply in the discrete dimension.  We handle this separately to avoid
    % potential problems with taking logarithms of functions with zeros.
    out = prodDiscrete(F);
end

end

function f = prodDiscrete(f)
%PRODDISCRETE   Compute product along discrete dimension.

% TODO:  We can probably do this more quickly by having the techs implement
% their own PROD and calling that in the array-valued case.  For now, we'll
% just convert to a quasimatrix and do it the slow way.
f = cheb2quasi(f);

p = f(1);
for k = 2:numel(f)
    p = p.*f(k);
end
f = p;

end
