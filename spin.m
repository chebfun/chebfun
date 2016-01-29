function [uout, tout] = spin(varargin)
%SPIN  Solve a time-dependent PDE in 1D with periodicity in space, using a 
%Fourier spectral method and an exponential integrator time-stepping scheme.
%
%   UOUT = SPIN(PDECHAR) solves the PDE specified by the string PDECHAR, and 
%   plots a movie of the solution as it computes it. The space and time 
%   intervals and the initial condition are chosen to produce beautiful movies. 
%   Strings available include 'AC' for Allen-Cahn equation, 'KS' for 
%   Kuramoto-Sivashinsky equation and 'KdV' for Korteweg-de Vries equation. 
%   Many other PDEs are available, see Remark 1 and Examples 1-10. The output 
%   UOUT is a CHEBFUN corresponding to the solution at the final time 
%   (a CHEBMATRIX for systems of equations, each row representing one variable). 
%
%   UOUT = SPIN(PDECHAR, TSPAN) solves the PDE from TPSAN(1) to TSPAN(END)
%   where TSPAN=[0 T1 T2 ... TF] is a vector of time chunks. The output UOUT is 
%   a CHEBMATRIX, each row corresponding to one variable and each column to one 
%   time chunk (unless TSPAN=[0 TF] and there is only one variable, in which 
%   case the output is a CHEBFUN at TF).
%
%   UOUT = SPIN(PDECHAR, TSPAN, U0) solves the PDE specified by the string 
%   PDECHAR on TSPAN x U0.DOMAIN, with initial condition a CHEBFUN U0 (one 
%   variable) or a CHEBMATRIX U0 (systems). See Example 11.
%
%   UOUT = SPIN(S, TSPAN, U0) solves the PDE specified by the spinop S and plots
%   a movie of the solution as it computes it. See HELP/SPINOP and Example 12.
%
%   UOUT = SPIN(..., PREF) allows one to use the preferences specified by the 
%   SPINPREF object PREF. See HELP/SPINPREF and Example 13.
% 
%   [UOUT, TOUT] = SPIN(...) also returns the times chunks TOUT at which UOUT
%   was computed.
%
% Remark 1: Available (case-insensitive) strings PDECHAR are
%
%    - 'AC' for Allen-Cahn equation, 
%    - 'Burg' for viscous Burgers equation,
%    - 'BZ' for Belousov-Zhabotinsky equation,
%    - 'CH' for Cahn-Hilliard equation,
%    - 'GS' for Gray-Scott equations,  
%    - 'KdV' for Korteweg-de Vries equation,
%    - 'KS' for Kuramoto-Sivashinsky equation, 
%    - 'Niko' for Nikolaevskiy equation, 
%    - 'NLS' for the focusing nonlinear Schroedinger equation,
%    - 'OK' for the Ohta-Kawasaki equation.
%
% Example 1: Allen-Cahn equation (metastable solutions)
%
%       u = spin('AC');
%
%    solves the Allen-Cahn equation 
%
%       u_t = 5e-3*u_xx - u^3 + u,
%
%    on [0 2*pi] from t=0 to t=300, with initial condition 
%
%       u0(x) = tanh(2*sin(x)) + 3*exp(-27*(x-4.2)^2) 
%               - 3*exp(-23.5*(x-pi/2)^2) + 3*exp(-38*(x-5.4)^2).
%
% Example 2: Viscous Burgers equation (shock formation and dissipation)
%
%       u = spin('Burg');
%
%    solves the viscous Burgers equation 
%
%           u_t = 1e-3*u_xx - u*u_x,
%
%   on [-1 1] from t=0 to t=20, with initial condition 
%
%       u0(x) = (1-x^2)*exp(-30*(x+1/2)^2.
%
% Example 3: Belousov-Zhabotinsky (reaction-diffusion with three species)
%
%       u = spin('BZ');
%
%    solves the Belousov-Zhabotinsky equation 
%
%       u_t = 1e-5*diff(u,2) + u + v - u*v - u^2,
%       v_t = 2e-5*diff(v,2) + w - v - u*v,
%       w_t = 1e-5*diff(w,2) + u - w
%
%    on [-1 1] from t=0 to t=30, with initial condition 
%
%       u0(x) = exp(-100*(x+.5)^2),
%       v0(x) = exp(-100*(x)^2),
%       w0(x) = exp(-100*(x-.5)^2).
%
% Example 4: Cahn-Hilliard equation (metastable solutions)
%
%       u = spin('CH');
%
%    solves the Cahn-Hilliard equation 
%
%       u_t = 1e-2*(-u_xx - 1e-3*u_xxxx + (u^3)_xx),
%   
%    on [-1 1] from t=0 to t=70, with initial condition 
%
%       u0(x) = (sin(4*pi*x))^5-sin(pi*x).
%
% Example 5: Gray-Scott equations (pulse splitting)
%
%       u = spin('GS');
%
%    solves the Gray-Scott equations 
%
%       u_t = diff(u,2) + 2e-2*(1-u) - u*v^2,
%       v_t = 1e-2*diff(u,2) - 8.62e-2*v + u*v^2,
%
%    on [-50 50] from t=0 to t=15000, with initial condition 
%
%       u0(x) = 1 - 1/2*sin(pi*(x-L)/(2*L))^100,
%       v0(x) = 1/4*sin(pi*(x-L)/(2*L))^100,
%           with L=50.
%
% Example 6: Korteweg-de Vries equation (two-soliton solution)
%
%       u = spin('KdV');
%
%    solves the Korteweg-de Vries equation 
%
%       u_t = -u*u_x - u_xxx,
%
%    on [-pi pi] from t=0 to t=0.03015, with initial condition 
%
%       u0(x) = 3*A^2*sech(.5*A*(x+2))^2 + 3*B^2*sech(.5*B*(x+1))^2,
%           with A=25 and B=16.
%
% Example 7: Kuramoto-Sivashinsky (chaotic attractor)
%
%       u = spin('KS');
%
%    solves the Kuramoto-Sivashinsky equation 
%
%       u_t = -u*u_x - u_xx - u_xxxx,
%
%    on [0 32*pi] from t=0 to t=300, with intial condition 
%
%       u0(x) = cos(x/16)*(1 + sin(x/16)).
%
% Example 8: Nikolaevskiy equation (chaotic attractor)
%
%       u = spin('Niko');
%
%    solves the Nikolaevskiy equation 
%
%       u_t = -u*u_x + 1e-1*u_xx + u_xxxx + u_xxxxxx,
%
%    on [0 32*pi] from t=0 to t=300, with intial condition 
%
%       u0(x) = cos(x/16)*(1 + sin(x/16)).
%
% Example 9: Nonlinear Schroedinger equation (breather solution)
%
%       u = spin('NLS');
%
%    solves the focusing Nonlinear Schroedinger equation 
%
%       u_t = i*u_xx + i*|u|^2*u,
%
%    on [-pi pi] from t=0 to t=20, with initial condition 
%
%       u0(x) = 2*B^2/(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*x)) - 1)*A,
%           with A=2 and B=1.
%
% Example 10: Ohta-Kawasaki equation (pattern formation)
%
%       u = spin('OK');
%
%    solves the Ohta-Kawasaki equation 
%
%       u_t = -u_xx - 1e-2*u_xxxx - 4(u - sum(u)) + (u^3)_xx,
%
%    on [0 2*pi] from t=0 to t=4, with initial condition 
%
%       u0(x) = cos(x)/2.
%
% Example 11: KdV with different time interval and intial condition
%
%       tspan = [0 1e-2];
%       u0 = chebfun(@(x) 1000*sin(x), [-pi pi]);
%       u = spin('KdV', tspan, u0);
%
%    solves the KdV equation on [-pi pi] from t=0 to t=1e-2, with initial
%    condition 
%
%       u0(x) = 1000*sin(x).
%
% Example 12: PDE specified by a SPINOP
%
%       L = @(u) -diff(u,2)-diff(u,4);
%       N = @(u) -.5*diff(u.^2);
%       S = spinop(L, N, [0 32*pi]);
%       u0 = chebfun(cos(x/16).*(1 + sin(x/16)), [0 32*pi], 'trig');
%       tspan = [0 300];
%       u = spin(S, u0, tspan);
%
%   is equivalent to u = spin('KS').
%
% Example 13: Using preferences
%
%       pref = spinpref('dt', 1e-5, 'N', 256, 'plot', 'waterfall');
%       u = spin('KdV', pref);
%
%   solves the KdV equation using a time-step dt=1e-5, N=256 grid points and 
%   produces a WATERFALL plot as opposed to a movie.
%
% See also SPINOP, SPINPREF, SPINSCHEME, SPIN2, SPIN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% SPIN is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end