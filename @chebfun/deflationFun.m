function out = deflationFun(Nu, u, r, p, alp, type)
% DEFLATIONFUN    Wrapper for CHEBMATRIX/DEFLATIONFUN
%
% See also chebmatrix.deflationFun

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the CHEBMATRIX method, which is what CHEBOP uses
out = deflationFun(Nu, u, chebmatrix(r), p, alp, type);

end