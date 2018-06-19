function f = randnfundisk(lambda)
%RANDNFUNDISK   Smooth random function on the unit disk
%   F = RANDNFUNDISK(LAMBDA) returns a DISKFUN of maximum
%   frequency about 2pi/LAMBDA and standard normal distribution N(0,1)
%   at each point.  F is obtained by restricting output of RANDNFUN2
%   to the unit disk.
%
%   RANDNFUNDISK() uses the default value LAMBDA = 1.
%
% Examples:
%
%   f = randnfundisk(0.5); plot(f)
%   colormap([0 0 0; 1 0 0]), caxis(norm(caxis,inf)*[-1 1]), axis off
%
% See also RANDNFUN, RANDNFUN2, RANDNFUNSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if nargin == 0
    lambda = 1;
end

if lambda > 1
    fsquare = randnfun2(lambda,1.25*[-1 1 -1 1]);         % correct in principle
else
    fsquare = randnfun2(lambda,1.25*[-1 1 -1 1],'trig');  % much faster
end

% Sample fsquare over a polar tensor product grid that is dense enough to
% resolve fsquare to machine precision.
[m,n] = length(fsquare);
m = max(m,n);

% Not sure about this multipliers:  Idea is that Chebyshev expansions
% require roughly pi/2 times the resolution as a Fourier expansion. The
% sqrt(2) comes from the maximum length of a Cartesian square grid cell.
% The 1.25 comes from the fact that the fsquare is over a grid that is 1.25
% times bigger than a grid from [-1,1] 
% (this number seems especially suspect).
n = round(1.25*m*sqrt(2)*pi/2);
n = n + 1-mod(n,2);
% Also this multiplier is probably not right for the Fourier expansion.
m = round(m*1.25*pi*sqrt(2)); m = m + mod(m,2);

% Create the tensor product grid
r = chebpts(n); r = r((n+1)/2:end); % r should be restricted to [0,1].
[tt,rr] = meshgrid(trigpts(m,[-pi pi]),r);

% Sample the random function over the square and create the diskfun.
F = fsquare(rr.*cos(tt),rr.*sin(tt));
f = diskfun(F);

end
