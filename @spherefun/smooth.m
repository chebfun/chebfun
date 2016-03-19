function g = smooth(f, dt)
%SMOOTH    Gaussian filtering on the sphere.  
%   G = SMOOTH(F), applies a low-pass filter to F. This is based on 
%   Gaussian filtering.  
% 
%   G = SMOOTH(F, dt), applies a low-pass filter to F with parameter dt.
%   The default to dt=1e-4, which is equivalent to running the Heat
%   equation for time dt with the initial condition being F. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 ) 
    dt = 1e-4; 
end
K = sqrt(1/dt)*1i;

% Find the length of f.
[n, m] = length(f); 

% Smoothing consists of solving the heat equation on the sphere, which is
% equivalent to applying a low-pass Gaussian filter.  To discretize the
% equation we use backward Euler:
% u^{n+1} = u^{n} + dt*L*u^{n+1}.
% which results in a Helmholtz equation with "wave-number" K = i*sqrt(1/dt)
%            L*u^{n+1) + K*u^{n+1} = -1/dt*u^{n}
% We only do one step of backward Euler to smooth the solution, with 
% u^{0} = f.
g = spherefun.Helmholtz(-1/dt*f, K, m, n);

end 