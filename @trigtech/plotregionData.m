function data = plotregionData(u, userEps, ignored)
%PLOTREGIONDATA   Useful data for plotting Chebfun strip.
%   PLOTREGIONDATA(U) returns a struct containing data that can be used for
%   plotting the Chebfun strip of U, which is symmetric about the real axis with 
%   half width log(1/EPS)/(pi*N).
%
%       boundary: boundary of the Chebfun ellipse in complex variable.
%       xlim: xlim for the plot.
%       ylim: ylim for the plot.
%       auxiliary: boundary defines a single period (plotted in dotted line). 
%
%   PLOTREGIONDATA(U, EPS) allows a user-specified EPS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    userEps = eps;
end

% Simplify first, to make sure length matches the epslevel:
u = simplify(u);

% Compute the intercept of the strip on the imaginary axis
M = length(u);
if ( mod(M,2) == 1 )
    N = (M-1)/2;
else
    N = M/2-1;
end
a = 1i*log(1/userEps)/(pi*N);
data.boundary = [-1 + a; 1 + a; nan; -1 - a; 1 - a];
data.xlim = [-1.1 1.1];
data.ylim = [-1.1*imag(a) 1.1*imag(a)];
data.auxiliary = [-1 + a; -1 - a; nan; 1 + a; 1 - a];

end
