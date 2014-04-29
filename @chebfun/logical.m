function f = logical(f)
%LOGICAL   CHEBFUN logical.
%   LOGICAL(F) returns a CHEBFUN which evaluates to one at all points where F is
%   non-zero and zero otherwise.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Add breaks at the roots (since these will take the value 0 in the output).
f = addBreaksAtRoots(f);

% Loop over the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = logical(f.funs{k});
end

% pointValues:
tol = vscale(f)*epslevel(f);
f.pointValues = abs(f.pointValues) > tol;

end
