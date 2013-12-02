function f = makeChebfun(u, dom)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
funs = cell(numel(u), 1);
for k = 1:numel(u)
    ct = chebtech2({[], flipud(u{k})});
    funs{k} = bndfun(ct, dom(k:k+1));
end
f = chebfun(funs);
u = chebmatrix({f});
end
