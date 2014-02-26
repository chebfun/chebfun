function f = toFunction(disc, values)
%TOFUNCTION Convert discrete values of an ULTRAS discretization to a CHEBFUN.

% ultraS represents a function using polynomial coefficients. This function
% converts those to a CHEBFUN. The input may be piecewise smooth, as indicated
% by the dimension property of the discretization.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

dom = disc.domain;          % Domain we're working on, including breakpoints
v = mat2cell(full(values), disc.dimension);  % break into smooth pieces
funs = cell(numel(v), 1);  % one fun per piece
for k = 1:numel(v)
    % Construct CHEBTECH2 objects from each pice
    ct = chebtech2({[], flipud(v{k})});
    % Assign each piece to a subinterval with a BNDFUN
    funs{k} = bndfun(ct, dom(k: k + 1));
end
% Conver the FUNS cell-array to a CHEBFUN.
f = chebfun(funs);
end
