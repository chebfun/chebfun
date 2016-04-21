function f = cumsum(f, dim)
%CUMSUM   Indefinite integral of a SPHEREFUN.
% 
%  F = CUMSUM(F) or F = CUMSUM(F, DIM) is not defined on a spherefun.
%
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a SPHEREFUN user.
% 
% See also SPHEREFUN/CUMSUM2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SPHEREFUN:CUMSUM:fail',...
      'The indefinite integral is not defined for a spherefun.');
  
end