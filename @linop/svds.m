function [U, S, V] = svds(L, k, sigma)
%SVDS  Find some singular values and vectors of a compact LINOP.
%   SVDS of a LINOP is currently not supported.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:LINOP:svds:noSupport', ...
    'LINOP/SVDS() is not currently supported.');

end
