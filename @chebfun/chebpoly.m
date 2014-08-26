function out = chebpoly(varargin)
%CHEBPOLY   Chebyshev polynomial coefficients.
%   CHEBPOLY(F) is deprecated. Please use CHEBCOEFFS().
%
% See also PLOTCOEFFS.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:CHEBFUN:chebpoly:deprecated', ...
    'CHEBPOLY is deprecated. Please use CHEBCOEFFS instead.');
warning('off', 'CHEBFUN:CHEBFUN:chebpoly:deprecated');

out = flipud(chebcoeffs(varargin{:})).';

end
