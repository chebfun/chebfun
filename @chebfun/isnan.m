function out = isnan(F)
%ISNAN   Test if a CHEBFUN is NaN.
%   ISNAN(F) returns TRUE if F has any NaN values and FALSE otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = false;

% Empty CHEBFUNs are not NaN.
if ( isempty(F) )
    return
end

out = zeros(1, numel(F));
for k = 1:numel(F)
    out(k) = columnIsnan(F);
end
out = any(out);

end

function out = columnIsnan(f)

out = false;

% Check for NaN pointValues:
if ( any(isnan(f.pointValues(:))) )
    out = true;
    return
end

% Check for NaN FUNs.
for k = 1:numel(f.funs)
    if ( isnan(f.funs{k}) )
        out = true;
        return
    end
end

end
