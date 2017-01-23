function varargout = cumprod(varargin)
%CUMPROD  Indefinite product integral of a DISKFUN. 
%   This is not allowed and returns an error.  This function exists so that
%   the error message is meaningful to a DISKFUN user.
%
% See also DISKFUN/CUMSUM, DISKFUN/SUM, DISKFUN/PROD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:DISKFUN:CUMPROD:notSupported',...
        'Indefinite product integral of a DISKFUN not supported.');

end
