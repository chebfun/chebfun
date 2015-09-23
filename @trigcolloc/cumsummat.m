function Q = cumsummat(N)
%CUMSUMMAT   Trigonometric Fourier integration matrix.
%   Q = CUMSUMMAT(N) is the matrix that maps function values at N equi-spaced
%   points to values of the integral of the interpolating trigonometric
%   polynomial at those points.

% [TODO]: Add support.
error('CHEBFUN:TRIGCOLLOC:cumsummat:notSupported', ...
    ['Indefinite integration is currently not supported for ' ...
    'TRIGCOLLOC discretization.\nPlease consider using ' ....
    'CHEBCOLLOC2 or ULTRAS discretization.']);   

end
