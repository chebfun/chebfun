function f = restrict(f, dom)
% RESTRICT  Restrict the domain of a SPHEREFUN.
%
% F = RESTRICT(F, DOM) is not supported for a spherefun F. A spherefun
% cannot be restricted. 

% Copyright 2016 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

error('CHEBFUN:SPHEREFUN:RESTRICT:fail',...
       'The domain of a spherefun cannot be restricted.')

end