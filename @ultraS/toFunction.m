function f = toFunction(disc, values)
%TOFUNCTION Convert discrete values to chebfun.

% ultraS represents a function using polynomial coefficients. This function
% converts those to a chebfun. The input may be piecewise smooth, as
% indicated by the dimension property of the discretization.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
dom = disc.domain;
v = mat2cell(full(values), disc.dimension);  % break into smooth pieces
funs = cell(numel(v), 1);  % one fun per piece
for k = 1:numel(v)
    ct = chebtech2({[], flipud(v{k})});
    funs{k} = bndfun(ct, dom(k: k + 1));
end
f = chebfun(funs);
end
