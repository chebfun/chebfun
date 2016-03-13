function g = smooth( f, dt )
% SMOOTH    Gaussian filtering on the sphere. 
% 
%  G = SMOOTH( F ), applies a low-pass filter to F. This is based on
%  Gaussian filtering.  
% 
%  G = SMOOTH( F, dt ), applies a low-pass filter to F with parameter dt.
%  The default to dt=1e-4, which is equivalent to running the Heat
%  equation for time dt with the initial condition being F. 

if ( nargin < 2 ) 
    dt = 1e-4; 
end
K = sqrt(1/dt)*1i;

% Find the length of f.
[m, n] = length( f ); 

% Solve the Helmholtz equation on the sphere to apply the Gaussian filter. 
% Since, g is expected to smoother than f, and mxn discretization should
% suffice. 
g = spherefun.Helmholtz( -1/dt*f, K, m, n);

end 