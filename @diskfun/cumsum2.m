function varargout = cumsum2(varargin)
%CUMSUM2   Double indefinite integral of a DISKFUN.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a DISKFUN user.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:DISKFUN:CUMSUM2:fail',...
      'The indefinite integral is not supported for a diskfun.');
end
