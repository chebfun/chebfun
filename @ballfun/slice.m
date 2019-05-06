function varargout = slice(varargin)
% SLICE(f, varargin) is not supported for ballfuns. 

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

error('CHEBFUN:BALLFUN:SLICE:fail',...
       'Slice is not supported for ballfuns. Use plot(f, ''slices'')')

end