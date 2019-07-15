function f = restrict(f, dom)
% RESTRICT  Restrict the domain of a DISKFUN.
%
% F = RESTRICT(F, DOM) is not supported for a diskfun F. A diskfun
% cannot be restricted. 

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

error('CHEBFUN:DISKFUN:RESTRICT:fail',...
       'The domain of a diskfun cannot be restricted.')

end