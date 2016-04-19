function varargout = quiver3(varargin)
%QUIVER3  3-D quiver plot of a SPHEREFUNV at data mapped by a SPHEREFUN.
%
%   This is not allowed and returns an error.  This function exists so that
%   the error message is meaningful to a SPHEREFUN user.
%
% See also SPHEREFUNV/QUIVER3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SPHEREFUN:QUIVER3:notSupported',...
        'quiver3 is not yet supported for SPHEREFUN.');

end
