function fout = cos(fin)
%COS Cosine of an ADCHEBFUN2.
%
% See also COSH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

fout = fin;
fout.chebfun2 = cos(fin.chebfun2);

fout.der = (-sin(fin.chebfun2))*fin.der;

end