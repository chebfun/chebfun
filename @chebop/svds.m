function [U, S, V] = svds(L, k, sigma)
%SVDS  Find some singular values and vectors of a compact linear CHEBOP.
%   SVDS of a CHEBOP is currently not supported.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CHEBOP:svds:nosupport', ...
    'CHEBOP/SVDS() is not currently supported.');

end
