function out = isperiodic(f)
%ISPERIODIC   Test if a FUN object is built upon TRIGTECH.
%   out = ISPERIODIC(F) returns logical true if its ONEFUN is a TRIGTECH
%   and false otherwise.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = isa(f.onefun, 'trigtech');

end
