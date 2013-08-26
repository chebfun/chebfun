function out = isempty(f)
%ISEMPTY   Test for empty chebfun.
%   ISEMPTY(F) returns logical true if F is an empty CHEBFUN and false
%   otherwise.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( (numel(f.funs) > 0) && ~isempty(f.funs(1)) )
    out = false;
else
    out = true;
end

end
