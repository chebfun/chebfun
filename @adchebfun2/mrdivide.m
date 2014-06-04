function F = mrdivide(f, g)
% /   Right scalar divide for ADCHEBFUN2.
% 
% F/C divides the chebfun2 F by a scalar C.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

F = f;
if (isa(g,'double'))
    F = F * ( 1 / g );
else
    error('ADCHEBFUN2:mrdivide:ADchebfun2','Did you mean ./ ?');
end

end