function f = convertToCorrectTech(f, newtech)
%CONVERTTOCORRECTTECH   Convert a CHEBFUN to another TECH.
%   CONVERTTOCORRECTTECH(F, NEWTECH) converts the CHEBFUN F to the TECH
%   NEWTECH.
%
% See also CHEBMATRIX/CONVERTTOCORRECTTECH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty or numeric case. 
if ( isempty(f) || isnumeric(f) )
    return
end

% If no TECH specified, do nothing.
if ( nargin == 1 )
    return
end

% Convert if necessary.
if ( ~isequal(get(f.funs{1}, 'tech'), newtech) )
    f = chebfun(f, f.domain, 'tech', newtech);
else
    return
end

end
