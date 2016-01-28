function [uout, tout] = spin3(varargin)
%SPIN3  Solve a time-dependent PDE in 3D with periodicity in space, using a 
%Fourier spectral method and an exponential integrator time-stepping scheme.
%
%   UOUT = SPIN3(PDECHAR) solves the PDE specified by the string PDECHAR, and 
%   plots a movie of the solution as it computes it. The space and time 
%   intervals, and the initial condition are chosen to produce beautiful movies. 
%   Strings available include include 'GL3' for Ginzburg-Landau equation and 
%   'GS3' for Gray-Scott equations. Many other PDEs are available, see Remark 1 
%   and Examples 1-4. The output UOUT is a CHEBFUN2 corresponding to the 
%   solution at the final time (a CHEBMATRIX for systems of equations, each row 
%   representing one variable). 
%
%   UOUT = SPIN3(PDECHAR, TSPAN) solves the PDE from TPSAN(1) to TSPAN(END)
%   where TSPAN=[0 T1 T2 ... TF] is a vector of time chunks. The output UOUT is 
%   a CHEBMATRIX, each row corresponding to one variable and each column to one 
%   time chunk (unless TSPAN=[0 TF] and there is only one variable, in that case 
%   the ouput is a CHEBFUN3 at TF).
%
%   UOUT = SPIN3(PDECHAR, TSPAN, U0) solves the PDE specified by the string 
%   PDECHAR on TSPAN x U0.DOMAIN, with initial condition a CHEBFUN3 U0 (one 
%   variable) or a CHEBMATRIX U0 (systems). 
%
%   UOUT = SPIN3(S, TSPAN, U0) solves the PDE specified by the spinop S and plots
%   a movie of the solution as it computes it. See HELP/SPINOP3.
%
%   UOUT = SPIN3(..., PREF) allows one to use the preferences specified by the 
%   SPINPREF object PREF. See HELP/SPINPREF3.
% 
%   [UOUT, TOUT] = SPIN3(...) also returns the times chunks TOUT at which UOUT
%   was computed.
%
% Remark 1: Available strings PDECHAR are
%
%    - 'GL3' for Ginzburg-Landau equation 
%
%           u_t = laplacian(u) + u - (1+1.3i)*u*|u|^2,
%
%    - 'GS3' for Gray-Scott equations
%
%           u_t = 2e-5*laplacian(u) + 3.5e-2*(1-u)*u - u*v^2,
%           v_t = 1e-5*laplacian(v) - 9.5e-2*v + u*v^2,
%           
%    - 'Schnak3' for Schnakenberg equations
%
%           u_t = laplacian(u) + .1 - u + u^2*v,
%           v_t = 10*laplacian(v) + .9 - u^2*v,
%
%    - and 'SH3' for Swift-Hohenberg equation
%
%           u_t = -2*laplacian(u) - biharmonic(u) - .9*u + u^2 - u^3.
%
% Example 1: Ginzburg-Landau equation (spiral waves)
%
%       u = spin('GL3');
%
%    solves the Ginzburg-Landau equation on [0 100]^3 from t=0 to t=200, with a
%    random initial condition.
%
% Example 2: Gray-Scott equations (pattern formation)
%
%       u = spin('GS3');
%
%    solves the Gray-Scott equations on [0 .75]^3 from t=0 to t=1600, with 
%    initial condition 
%
%       u0(x,y,z) = 1 - exp(-150*((x-G/2).^2 + (y-G/2).^2 + (z-G/2).^2)),
%       v0(x,y,z) = exp(-150*((x-G/2).^2 + 2*(y-G/2).^2 + (z-G/2).^2)),
%           with G=.75.
%
% Example 3: Schnakenberg equations (pattern formation)
%
%       u = spin('Schnak3');
%
%    solves the Schnakenberg equations on [0 20]^3 from t=0 to t=400, with 
%    initial condition 
%
%       u0(x,y,z) = 1 - exp(-10*((x-G/2).^2 + (y-G/2).^2 + (z-G/2).^2)),
%       v0(x,y,z) = exp(-10*((x-G/2).^2 + 2*(y-G/2).^2 + (z-G/2).^2)),
%           with G=20.
%
% Example 4: Swift-Hohenberg equation (pattern formation)
%
%       u = spin('SH3');
%
%    solves the Swift-Hohenberg equation on [0 50]^3 from t=0 to t=200, with 
%    initial condition
%
%       u0(x,y,z) = 1/4*(sin(pi*x/10) + sin(pi*y/10) + sin(pi*z/10)
%                       + sin(pi*x/2).*sin(pi*y/2) + sin(pi*x/2).*sin(pi*z/2)
%                       + sin(pi*z/2).*sin(pi*y/2)).
%
% See also SPINOP3, SPINPREF3, SPINSCHEME, SPIN, SPIN2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% SPIN3 is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end