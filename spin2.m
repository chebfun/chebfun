function [uout, tout] = spin2(varargin)
%SPIN2  Solve a time-dependent PDE in 2D with periodicity in space, using a 
%Fourier spectral method and an exponential integrator time-stepping scheme.
%
%   UOUT = SPIN2(PDECHAR) solves the PDE specified by the string PDECHAR, and 
%   plots a movie of the solution as it computes it. The space and time 
%   intervals and the initial condition are chosen to produce beautiful movies. 
%   Strings available include include 'GL2' for Ginzburg-Landau equation and 
%   'GS2' for Gray-Scott equations. Many other PDEs are available, see Remark 1 
%   and Examples 1-4. The output UOUT is a CHEBFUN2 corresponding to the 
%   solution at the final time (a CHEBMATRIX for systems of equations, each row 
%   representing one variable). 
%
%   UOUT = SPIN2(PDECHAR, TSPAN) solves the PDE from TPSAN(1) to TSPAN(END)
%   where TSPAN=[0 T1 T2 ... TF] is a vector of time chunks. The output UOUT is 
%   a CHEBMATRIX, each row corresponding to one variable and each column to one 
%   time chunk (unless TSPAN=[0 TF] and there is only one variable, in which 
%   case the output is a CHEBFUN2 at TF).
%
%   UOUT = SPIN2(PDECHAR, TSPAN, U0) solves the PDE specified by the string 
%   PDECHAR on TSPAN x U0.DOMAIN, with initial condition a CHEBFUN2 U0 (one 
%   variable) or a CHEBMATRIX U0 (systems). 
%
%   UOUT = SPIN2(S, TSPAN, U0) solves the PDE specified by the spinop S and plots
%   a movie of the solution as it computes it. See HELP/SPINOP2.
%
%   UOUT = SPIN2(..., PREF) allows one to use the preferences specified by the 
%   SPINPREF object PREF. See HELP/SPINPREF2.
% 
%   [UOUT, TOUT] = SPIN2(...) also returns the times chunks TOUT at which UOUT
%   was computed.
%
% Remark 1: Available (case-insensitive) strings PDECHAR are
%
%    - 'GL2' for Ginzburg-Landau equation, 
%    - 'GS2' for Gray-Scott equations,
%    - 'Schnak2' for Schnakenberg equations,
%    - 'SH2' for Swift-Hohenberg equation.
%
% Example 1: Ginzburg-Landau equation (spiral waves)
%
%       u = spin2('GL2');
%
%    solves the Ginzburg-Landau equation 
%
%        u_t = laplacian(u) + u - (1+1.3i)*u*|u|^2,
%
%    on [0 200]^2 from t=0 to t=80, with a random initial condition.
%
% Example 2: Gray-Scott equations (fingerprints patterns)
%
%       u = spin2('GS2');
%
%    solves the Gray-Scott equations 
%
%       u_t = 2e-5*laplacian(u) + 3.5e-2*(1-u)*u - u*v^2,
%       v_t = 1e-5*laplacian(v) - 9.5e-2*v + u*v^2,
%
%    on [0 1.25]^2 from t=0 to t=3200, with initial condition 
%
%       u0(x,y) = 1 - exp(-150*((x-G/2)^2 + (y-G/2)^2)),
%       v0(x,y) = exp(-150*((x-G/2)^2 + 2*(y-G/2)^2)),
%           with G=1.25.
%
% Example 3: Schnakenberg equations (pattern formation)
%
%       u = spin2('Schnak2');
%
%    solves the Schnakenberg equations 
%
%       u_t = laplacian(u) + .1 - u + u^2*v,
%       v_t = 10*laplacian(v) + .9 - u^2*v,
%
%    on [0 50]^2 from t=0 to t=200, with initial condition 
%
%       u0(x,y) = 1 - exp(-10*((x-G/2)^2 + (y-G/2)^2))
%       v0(x,y) = exp(-10*((x-G/2)^2 + 2*(y-G/2)^2)),
%           with G=50.
%
% Example 4: Swift-Hohenberg equation (Rayleigh-Benard convection)
%
%       u = spin2('SH2');
%
%    solves the Swift-Hohenberg equation 
%
%       u_t = -2*laplacian(u) - biharmonic(u) - .9*u + u^2 - u^3,
%
%    on [0 50]^2 from t=0 to t=200, with a random initial condition.
%
% See also SPINOP2, SPINPREF2, SPINSCHEME, SPIN, SPIN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

for j = 1:nargin
    item =  varargin{j};
    if ( isa(item, 'spinoperator') == 1 )
        if ( isa(item, 'spinop') == 1 )
            error('CHEBFUN:SPIN', 'Use SPIN for PDEs in one space dimension.')
        elseif ( isa(item, 'spinop3') == 1 )
            error('CHEBFUN:SPIN', ['Use SPIN3 for PDEs in three space ', ...
                'dimensions.'])
        end
    elseif ( isa(item, 'char') == 1 )
        is1D = isempty(strfind(item, '2')) && isempty(strfind(item, '3'));
        is3D = ~isempty(strfind(item, '3'));
        if ( is1D == 1 )
            error('CHEBFUN:SPIN', 'Use SPIN for PDEs in one space dimension.')
        elseif ( is3D == 1 )
            error('CHEBFUN:SPIN', ['Use SPIN3 for PDEs in three space ', ...
                'dimensions.'])
        end
    end
end

% SPIN2 is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end