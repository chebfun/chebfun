function [U, S, V] = svds(L, k)
%SVDS  Find some singular values and vectors of a compact linear CHEBOP.
%   SVDS of a CHEBOP is currently not supported.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%error('CHEBFUN:CHEBOP:svds:nosupport', ...
%    'CHEBOP/SVDS() is not currently supported.');

% Linearize
[ K, ~, isLinear ] = linearize( L, [], [], 1 );
K

% check linearity
if ( any(isLinear == 0) )
    error('CHEBFUN:CHEBOP:svds:linear', ...
        'SVDS can only be applied to linear chebops.');
end

% construct adjoint
[ U, S, V ] = svds( K, k, getBCType(L) );

end
