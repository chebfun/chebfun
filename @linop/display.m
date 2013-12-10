function display(L)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

[m, n] = size(L);
fprintf('\n  %ix%i linear operator\n', m,n)

nc = size(L.constraint,1);
if ( nc == 1 )
    fprintf('\n    with 1 constraint/boundary condition\n')
elseif ( nc > 0 )
    fprintf('\n    with %i constraints/boundary conditions\n',nc)
end
fprintf('\n')

end
