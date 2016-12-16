function out = isPeriodicTech(f)
%ISPERIODICTECH   Test if a CHEBFUN object is built upon a TECH of 
%periodic functions.
%   out = ISPERIODICTECH(F) returns logical true if F has at least one FUN
%   which is made of a TECH of periodic functions and false otherwise.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% out = all(cellfun(@(f) isPeriodicTech(f), f.funs));

% Since all funs have the same technology, it is sufficient and faster to test
% only the first.
out = isPeriodicTech(f.funs{1});

end