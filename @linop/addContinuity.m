function L = addContinuity(L, varargin)
%ADDCONSTRAINT Append to linop constraints.
%   L = ADDCONTINUITY(L,FUN,VAL) adds a new constraint on the linop L. The
%   functional FUN when applied to a function will be required to be VAL.
%
% See also LINOP.ADDCONSTRAINT, LINOPCONSTRAINT. 

%  Copyright 2016 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% We don't really care about the content. But because the continuity property is
% modified, other linop calls will not apply automatically derived continuity
% conditions.

% Append the input constraint to the LINOPCONSTRAINT currently stored in the
% input linop.
L.continuity = append(L.continuity, varargin{:});

% We may have to update the domain, if a new breakpoint was introduced.
n1 = length(L.domain);
F = L.continuity.functional;
n2 = length(F.domain);
if ( n1 ~= n2 )
    L.domain = domain.merge(L.domain,F.domain);
end

end
