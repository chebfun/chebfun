function F = mrdivide(f,g)
% /   Right scalar divide
% 
% F/C divides the chebfun2 F by a scalar C.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

F = f;
if (isa(g,'double'))
    F.fun2 = F.fun2/g;
    F.scl = F.scl/g;
else
    error('CHEBFUN2:mrdivide:chebfun2chebfun2','Did you mean ./ ?');
end

end