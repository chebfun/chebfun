function data = plotregionData(u, userEps, numpts)
%PLOTREGIONDATA   Useful data for plotting Chebfun ellipse.
%   PLOTREGIONDATA(U) returns a struct containing data that can be used for
%   plotting the Chebfun ellipse of U, whose semi-minor and major axes summing
%   to rho = exp(abs(log(EPS))/length(U)), where EPS is the machine precision.
%   The struct DATA contains the following fields:
%
%       boundary: boundary of the Chebfun ellipse in complex variable.
%       xlim: xlim for the plot.
%       ylim: ylim for the plot.
%       auxiliary: empty. 
%
%   PLOTREGIONDATA(U, EPS) allows a user-specified EPS.
%
%   PLOTREGIONDATA(U, EPS, NUMPTS) uses N points to represent the boundary.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    numpts = 500;
end

if ( nargin < 2 )
    userEps = eps;
end

% Simplify first, to make sure length matches the epslevel:
u = simplify(u);

% The unit circle.
c = exp(2*pi*1i*linspace(0, 1, numpts).');
rho = exp(abs(log(userEps)) / length(u));
data.boundary = .5*(rho*c + 1./(rho*c));
a = (rho+1/rho);
b = (rho-1/rho);
data.xlim = [-1.1*a 1.1*a];
data.ylim = [-1.1*b 1.1*b];
data.auxiliary = [];

end
