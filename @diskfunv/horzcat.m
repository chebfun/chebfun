function varargout = horzcat( varargin )
%HORZCAT   Horizontal concatenation of CHEBFUN2V objects.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a CHEBFUN2 user.

% Copyright 2015 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

error('CHEBFUN:DISKFUNV:horzcat:notSupported', ...
    'Horizontal concatenation of DISKFUNV objects is not supported.')

end
