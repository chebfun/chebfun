function f = fliplr(f)
%FLIPUD  Flip/reverse a chebfun.
%   G = FLIPLR(F), where F is a row CHEBFUN, returns a chebfun G with the same
%   domain as F but reversed; that is, G(x) = F(a+b-x), where the domain is
%   [a,b].
%
%   FLIPLR(F), where F is an array-valued column chebfun, exchanges the order of
%   the columns of F.
%
% See also CHEBFUN/FLIPUD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

if ( ~f.isTransposed )

    % Flip the coluns of an array-valued chebfun:
    for k = 1:numel(f.funs)
        f.funs{k} = fliplr(f.funs{k});
    end

else

    % Transpose f and call FLIPUD():
    f = transpose(flipud(transpose(f)));

end