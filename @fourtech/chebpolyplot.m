function varargout = chebpolyplot(f, varargin)
%CHEBPOLYPLOT   Display Fourier coefficients graphically.
%
% See also COEFPLOT, PLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This is just a wrapper to make things work with chebfun.  Simply call
% coefplot to produce the desired result.

% Give an output if one was requested:
varargout = [];
if ( nargout > 0 )
    varargout{1} = coefplot(f,varargin{:});
else
    coefplot(f,varargin{:});
end

end
