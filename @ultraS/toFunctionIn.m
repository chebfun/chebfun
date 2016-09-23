function f = toFunctionIn(disc, coeffs)
%TOFUNCTION   Convert discrete values of an ULTRAS discretization to a CHEBFUN.
%   F = TOFUNCTIONIN(DISC, COEFFS) converts the coeffs returned by ULTRAS to a
%   CHEBFUN. The input may be piecewise smooth, as indicated by the dimension
%   property of the discretization.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

dom = disc.domain;         % Domain we're working on, including breakpoints
c = mat2cell(full(coeffs), disc.dimension); % Break into smooth pieces
funs = cell(numel(c), 1);  % One FUN per piece
for k = 1:numel(c)
    % Construct CHEBTECH2 objects from each piece:
    ct = chebtech2({[], c{k}});
    % Assign each piece to a subinterval with a BNDFUN:
    funs{k} = bndfun(ct, struct('domain', dom(k: k + 1)));
end
% Conver the FUNS cell-array to a CHEBFUN.
f = chebfun(funs);

end
