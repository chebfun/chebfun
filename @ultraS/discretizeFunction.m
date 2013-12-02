function f_coeffs = discretizeFunction(f, dim, dom)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
if ( nargin < 3 )
    dom = f.domain;
end
f_coeffs = [];
f = restrict(f, dom);
for k = 1:numel(dom)-1
    dimk = dim(k);
    tmp = flipud(get(f.funs{k}, 'coeffs'));
    n = length(tmp);
    % prolong/truncate.
    if ( n > dimk )
        tmp = tmp(1:dimk);
    else
        tmp = [tmp ; zeros(dimk - n, 1)];
    end
    f_coeffs = [f_coeffs ; tmp];
end
end
