function varargout = horzcat( varargin )
%HORZCAT   Horizontal concatenation of DISKFUNV objects.
%   This is not allowed and returns an error.  This function exists so that 
%   the error message is meaningful to a DISKFUNV user.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:DISKFUNV:horzcat:notSupported', ...
    'Horizontal concatenation of DISKFUNV objects is not supported.')

end