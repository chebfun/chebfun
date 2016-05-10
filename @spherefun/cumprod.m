function varargout = cumprod(varargin)
%CUMPROD  Indefinite product integral of a SPHEREFUN. 
%   This is not allowed and returns an error.  This function exists so that
%   the error message is meaningful to a SPHEREFUN user.
%
% See also SPHEREFUN/CUMSUM, SPHEREFUN/SUM, SPHEREFUN/PROD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SPHEREFUN:CUMPROD:notSupported',...
        'Indefinite product integral of a SPHEREFUN not supported.');

end
