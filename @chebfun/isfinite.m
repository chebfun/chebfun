function out = isfinite(F)
%ISFINITE   Test if a CHEBFUN is bounded.
%   ISFINITE(F) returns FALSE if F has any infinite values and TRUE otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = zeros(1, numel(F));
for k = 1:numel(F)
    out(k) = columnIsfinite(F(k));
end
out = any(out);

end

function out = columnIsfinite(f)

% Check for infinite breakpoint values:
if ( any(~isfinite(f.pointValues(:))) )
    out = false;
    return
end

% Check for any infinite FUNs.
for k = 1:numel(f.funs)
    if ( ~isfinite(f.funs{k}) )
        out = false;
        return
    end
end

out = true;

end
