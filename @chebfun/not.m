function f = not(f)
%~   CHEBFUN logical NOT.
%   NOT(F) returns a CHEBFUN which evaluates to zero at all points where F is
%   zero and one otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Add breaks at the roots (since these will take the value 1 in the output).
f = addBreaksAtRoots(f);

% Loop over the FUNs:
for k = 1:numel(f.funs)
    f.funs{k} = not(f.funs{k});
end

% pointValues:
tol = vscale(f)*epslevel(f);
f.pointValues = abs(f.pointValues) < tol;

end
