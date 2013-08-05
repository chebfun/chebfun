function out = isreal(f)
%ISREAL   True for real-valued CHEBFUN object.
%   ISREAL(F) returns logical true if F does not have an imaginary part and
%   false otherwise.
%
%   ~ISREAL(F) detects chebfuns that have an imaginary part even if it is all
%   zero.
%
% See also REAL, IMAG.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check to see is the impulses are real:
out = isreal(f.impulses);
if ( ~out )
    % Complex number found. Break:
    return
end

% Check to see if each of the funs are real:
for k = 1:numel(f.funs)
    out = isreal(f.funs{k});
    if ( ~out )
        % Complex number found. Break:
        return
    end
end

end
