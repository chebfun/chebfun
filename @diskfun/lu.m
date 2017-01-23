function varargout = lu(varargin)
%LU   LU factorization of a DISKFUN.
%
%   This is not allowed and returns an error.  This function exists so that
%   the error message is meaningful to a DISKFUN user.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:DISKFUN:LU:notSupported',...
        'LU factorization of a DISKFUN is not supported.');

end
