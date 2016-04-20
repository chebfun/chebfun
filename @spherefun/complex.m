function varargout = complex(varargin)
%COMPLEX Construct complex SPHEREFUN from real and imaginary parts.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a SPHEREFUN user.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SPHEREFUN:COMPLEX:notSupported',...
        'Complex-valued spherefuns are not supported.');
end
