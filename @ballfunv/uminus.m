function g = uminus(f)
%- Unary minus of a BALLFUNV.
%   -F returns the BALLFUNV negated componentwise. 
%
%   UMINUS(F) is called by the syntax -F. 
%
% See also UPLUS.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F = f.comp;
g = [ -F{1} ; -F{2} ; -F{3} ];
end
