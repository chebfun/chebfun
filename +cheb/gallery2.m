function [f,fun] = gallery2(N)
%GALLERY2   Gallery of 2-dimensional functions.
%   GALLERY2(N) returns an interesting 2D function as a CHEBFUN2. N must
%   be between 1 and 7.
%
%   [F,FA] = GALLERY2(N) also returns an anonymous function used to define
%   the CHEBFUN2.
%

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F{1} = @(x,y) 1./(1+100*(x.^2-y.^2).^2);
F{2} = @(x,y) 1./(1+1e3*((x.^2-.25).^2.*(y.^2-.25).^2));
F{3} = @(x,y) cos(10*(x.^2+y)).*sin(10*(x+y.^2));
F{4} = @(x,y) real(airy(5*(x+y.^2)).*airy(-5*(x.^2+y.^2)));
F{5} = @(x,y) tanh(10*x).*tanh(10*y)./tanh(10).^2+cos(5*x);
F{6} = @(x,y) 3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ...
         - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ...
         - 1/3*exp(-(x+1).^2 - y.^2);  % PEAKS function

domain{1} = [-1 1 -1 1];
domain{2} = [-1 1 -1 1];
domain{3} = [-1 1 -1 1];
domain{4} = [-1 1 -1 1];
domain{5} = [-1 1 -1 1];
domain{6} = 3*[-1 1 -1 1];

f = chebfun2(F{N}, domain{N});
fun = F{N};

end