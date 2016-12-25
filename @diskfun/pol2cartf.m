function g = pol2cartf(f,th,r)
%POL2CARTF     Wrapper for evaluating f(x, y) in polar coordinates. 
% 
% G = POL2CARTF( F, TH, R ) converts the function F(X, Y) to a function G so
% that G(TH, R) = F(R.*cos(TH), R.*sin(TH)). 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
 
g = feval(f, r.*cos(th), r.*sin(th));

end