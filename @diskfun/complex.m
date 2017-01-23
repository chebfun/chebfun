function varargout = complex(varargin)
% COMPLEX  Construct complex DISKFUN from real and imaginary parts.
%   This is not allowed and returns an error.  This function exists so that
%   the error message is meaningful to a DISKFUN user.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:DISKFUN:COMPLEX:notSupported',...
        'Complex-valued diskfuns are currently not supported.');
end
