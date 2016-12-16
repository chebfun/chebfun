function f = cumsum(f, dim)
%CUMSUM   Indefinite integral of a DISKFUN.
% 
%  F = CUMSUM(F) or F = CUMSUM(F, DIM) is not defined on a diskfun.
%
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a DISKFUN user.
% 
% See also DISKFUN/CUMSUM2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:DISKFUN:CUMSUM:fail',...
      'The indefinite integral is not defined for a diskfun.');
  
end