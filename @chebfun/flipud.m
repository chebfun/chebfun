function f = flipud(f)
%FLIPUD   Flip/reverse a CHEBFUN.
%   G = FLIPUD(F), where F is a column CHEBFUN, returns a CHEBFUN G with the
%   same domain as F but reversed; that is, G(x) = F(a+b-x), where the domain is
%   [a,b].
%
%   FLIPUD(F), where F is an array-valued row CHEBFUN, exchanges the order of
%   the rows of F.
%
% See also CHEBFUN/FLIPLR.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

if ( ~f.isTransposed )

    % Reverse and translate the breakpoints.
    newDomain = -fliplr(f.domain) + sum(f.domain([1, end]));
    % Reverse the order of the corresponding impulses:
    f.impulses = fliplr(f.impulses);

    % Reverse the order of funs:
    f.funs = f.funs(end:-1:1);
    % and the funs themselves.
    for k = 1:numel(f.funs)
        f.funs{k} = flipud(f.funs{k});
        f.funs{k} = map(f.funs{k}, newDomain(k:k+1));
    end

    % Apply the new domain to the chebfun:
    f.domain = newDomain;

else

    f = transpose(flipud(transpose(f)));

end
