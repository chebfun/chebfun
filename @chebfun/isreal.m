function out = isreal(F)
%ISREAL   True for real-valued CHEBFUN object.
%   ISREAL(F) returns logical true if F does not have an imaginary part and
%   false otherwise.
%
%   Unlike the built in MATLAB function, ~ISREAL(F) does not detect CHEBFUN
%   objects that have an all zero imaginary part.
%
% See also REAL, IMAG.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = zeros(1, numel(F));
for k = 1:numel(F)
    out(k) = columnIsreal(F(k));
end
out = all(out);

end

function out = columnIsreal(f)

% Check to see is the impulses are real:
out = isreal(f.pointValues);
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
