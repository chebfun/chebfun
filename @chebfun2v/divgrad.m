function f = divgrad(f)
%DIVGRAD   Laplacian of a CHEBFUN2V.
%   F = DIVGRAD(F) returns the Laplacian of a CHEBFUN2V i.e.,
%       divgrad(f) = f_xx + f_yy 
%
% Also see CHEBFUN2V/LAP.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

f = diff(f.xcheb, [2,0]) + diff(f.ycheb, [0,2]);
 
end