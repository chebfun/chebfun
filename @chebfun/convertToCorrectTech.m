function f = convertToCorrectTech(f, domain, newtech)
%CONVERTTOCORRECTTECH   Convert a CHEBFUN to another TECH.
%   CONVERTTOCORRECTTECH(F, NEWTECH) converts the CHEBFUN F to the TECH
%   NEWTECH.
%
% See also CHEBMATRIX/CONVERTTOCORRECTTECH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty case.
if ( isempty(f) )
    return
end

% If no TECH specified, do nothing.
if ( nargin == 1 )
    return
end

% Numeric.
if ( isnumeric(f) )
   f = chebfun(f, domain, 'tech', newtech);
   return
end

% Convert if necessary.
if ( ~isequal(get(f.funs{1}, 'tech'), newtech) )
    f = chebfun(f, domain, 'tech', newtech);
else
    return
end

end
