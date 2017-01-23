function F = removeDeltas(F)
%REMOVEDELTAS   Remove the deltas from a chebfun.
%   REMOVEDELTAS(F) returns a chebfun identical to F with all deltas removed.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

for k = 1:length(F.funs)
    F.funs{k} = removeDeltas(F.funs{k});
end

F.pointValues = chebfun.getValuesAtBreakpoints(F.funs);

end
