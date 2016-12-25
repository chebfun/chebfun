function varargout = horzcat( varargin )
%HORZCAT   Horizontal concatenation of CHEBFUN2V objects.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a CHEBFUN2 user.

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

error('CHEBFUN:CHEBFUN2V:horzcat:notSupported', ...
    'Horizontal concatenation of CHEBFUN2V objects is not supported.')

end
