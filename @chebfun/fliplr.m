function F = fliplr(F)
%FLIPLR   Flip/reverse a CHEBFUN.
%   G = FLIPLR(F), where F is a row CHEBFUN, returns a CHEBFUN G with the same
%   domain as F but reversed; that is, G(x) = F(a+b-x), where the domain is
%   [a,b].
%
%   FLIPLR(F), where F is an array-valued column CHEBFUN or a quasimatrix,
%   reverses the order of the columns of F.
%
% See also FLIPUD.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) )
    return
end

if ( numel(F) > 1 && ~F(1).isTransposed )
    idx = fliplr(1:numColumns(F));
    F = F(idx);
else
    for k = 1:numel(F)
        F(k) = columnFliplr(F(k));
    end
end

end

function f = columnFliplr(f)

if ( ~f.isTransposed )

    % Flip the columns of an array-valued chebfun:
    for k = 1:numel(f.funs)
        f.funs{k} = fliplr(f.funs{k});
    end

    % Flip the pointValues:
    f.pointValues = fliplr(f.pointValues);

else

    % Transpose f and call FLIPUD():
    f = flipud(f.').';

end

end
