function varargout = horzcat(varargin)
%HORZCAT   Horizontal concatenation of SPHEREFUNV objects.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a SPHEREFUN user.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('SPHEREFUN:SPHEREFUNV:horzcat:notSupported', ...
    'Horizontal concatenation of SPHEREFUNV objects is not supported.')

end
