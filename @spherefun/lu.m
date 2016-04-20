function varargout = lu(varargin)
%LU   LU factorization of a SPHEREFUN.
%
%   This is not allowed and returns an error.  This function exists so that
%   the error message is meaningful to a SPHEREFUN user.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SPHEREFUN:CUMPROD:notSupported',...
        'LU factorization of a SPHEREFUN is not supported.');

end
