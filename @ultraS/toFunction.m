function f = toFunction(disc, values)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
dom = disc.domain;
v = mat2cell(full(values), disc.dimension);
funs = cell(numel(v),1);
for k = 1:numel(v)
    ct = chebtech2({[], flipud(v{k})});
    funs{k} = bndfun(ct, dom(k:k+1));
end
f = chebfun(funs);
end
