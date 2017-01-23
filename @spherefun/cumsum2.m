function varargout = cumsum2(varargin)
%CUMSUM2   Double indefinite integral of a SPHEREFUN.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a SPHEREFUN user.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SPHEREFUN:CUMSUM2:fail',...
      'The indefinite integral is not supported for a spherefun.');
end
