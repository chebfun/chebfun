function F = minus( F, G )
% - MINUS. Minus of two DISKFUNV.  
%   F - G subtracts the DISKFUNV F from G componentwise.
%
%   minus(F, G) is called for the syntax f - g.
% 
% See also DISKFUNV/PLUS

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% F - G = F + (-G):
F = plus( F, uminus( G ) );

end