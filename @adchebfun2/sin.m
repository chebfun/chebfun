function f = sin( f ) 
% SIN   Sine of an ADCHEBFUN2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

fout = fin;
fout.chebfun2 = sin(fin.chebfun2);

fout.der = (cos(fin.chebfun2))*fin.der;

end