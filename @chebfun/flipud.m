function F = flipud(F)
%FLIPUD   Flip/reverse a CHEBFUN.
%   G = FLIPUD(F), where F is a column CHEBFUN, returns a CHEBFUN G with the
%   same domain as F but reversed; that is, G(x) = F(a+b-x), where the domain is
%   [a,b].
%
%   FLIPUD(F), where F is an array-valued row CHEBFUN or a quasimatrix, reverses
%   the order of the rows of F.
%
% See also FLIPLR.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) )
    return
end

if ( numel(F) > 1 && F(1).isTransposed )
    idx = fliplr(1:numColumns(F));
    F = F(idx);
else
    for k = 1:numel(F)
        F(k) = columnFlipud(F(k));
    end
end

end
   
function f = columnFlipud(f)

if ( ~f.isTransposed )

    % Reverse and translate the breakpoints.
    newDomain = -fliplr(f.domain) + sum(f.domain([1, end]));
    % Reverse the order of the corresponding pointValues:
    f.pointValues = flipud(f.pointValues);

    % Reverse the order of FUNs:
    f.funs = f.funs(end:-1:1);
    % and the FUNs themselves.
    for k = 1:numel(f.funs)
        f.funs{k} = flipud(f.funs{k});
        f.funs{k} = changeMap(f.funs{k}, newDomain(k:k+1));
    end

    % Apply the new domain to the chebfun:
    f.domain = newDomain;

else

    % Transpose f and call FLIPLR():
    f = fliplr(f.').';

end

end
