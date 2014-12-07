function f = changeTech(f, newtech)
%CHANGETECH   Convert a CHEBFUN to another TECH.
%   CHANGETECH(F, NEWTECH) converts the CHEBFUN F to the TECH NEWTECH.
%
% See also CHEBMATRIX/CHANGETECH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty case. 
if ( isempty(f) )
    return
end

% Convert if necessary.
if ( ~isequal(get(f.funs{1}, 'tech'), newtech) )
    f = chebfun(f, f.domain, 'tech', newtech);
else
    return
end

end