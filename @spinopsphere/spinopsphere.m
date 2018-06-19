classdef spinopsphere < spinoperator
%SPINOPSPHERE   Class for representing the spatial part of differential
%operators for time-dependent PDEs on the sphere.
%   SPINOPSPHERE is a class for representing the spatial part S of a 
%   time-dependent PDE of the form u_t = S(u) = Lu + N(u) on the sphere, where 
%   L is a linear operator and N is a nonlinear operator. 
%
%   S = SPINOPSPHERE(PDECHAR) creates a SPINOPSPHERE object S defined by the 
%   string PDECHAR. Strings available include 'AC' for Allen-Cahn equation, 'GL'
%   GL' for Ginzburg-Landau equation, 'GM' for Gierer-Meinhardt equations and 
%   'NLS' for the nonlinear Schroedinger equation.
%
%   S = SPINOPSPHERE(TSPAN) creates a SPINOPSPHERE object S on SPHERE x TSPAN. 
%   The other fields of a SPINOPSPHERE are its linear part S.LIN, its nonlienar 
%   part S.NONLIN and the initial condition S.INIT. The fields can be set via 
%   S.PROP = VALUE. See Remarks 1/2 and Example 1.
%
% Remark 1: The linear part has to be of the form 
%           
%               @(u) A*lap(u)
%
%           for some number A.
%
% Remark 2: The nonlinear part has to be of the form @(u) f(u), where f is a 
%           nonlinear nondifferential operator with constant coefficients.
%
% Example 1: To construct a SPINOPSPHERE corresponding to the GL equation on 
%            the sphere with TSPAN = [0 100] and initial condition 
%            u0(lam,th) = Y_8^2(lam,th), one can type
%
%            tspan = [0 100];
%            S = spinopsphere(tspan);
%            S.lin = @(u) 5e-4*lap(u);
%            S.nonlin = @(u) u - (1+1.5i)*u.*(abs(u).^2);
%            S.init = spherefun.sphharm(8, 2);
%
% See also SPINOPERATOR, SPINSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinopsphere(varargin)
            
            % Zero input argument:
            if ( nargin == 0 )
                return
                        
            % One input argument:
            elseif ( nargin == 1 )
                            
                S.domain = [-pi pi 0 pi];
                            
                if ( isa(varargin{1}, 'char') == 1 )
                    [L, N, tspan, u0] = parseInputs(varargin{1});
                    S.lin = L;
                    S.nonlin = N;
                    S.tspan = tspan;
                    S.init = u0;
                elseif ( isa(varargin{1}, 'double') == 1 )
                    S.tspan = varargin{1};
                else
                    error('SPINOPSPHERE:constructor', ['When constructing ', ...
                        'a SPINOPSPHERE with one input argument, this ', ...
                        'argument should be a STRING or a DOUBLE.'])
                end
         
            % More than one input argument:    
            else
                error('SPINOPSPHERE:constructor', 'Too many input arguments.')
            end
            
        end
        
    end

end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% METHODS IMPLEMENTED IN THIS FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    function [L, N, tspan, u0] = parseInputs(pdechar)
        %PARSEINPUTS   Parse inputs when using a STRING.
        %   [L, N, TSPAN, U0] = PARSEINPUTS(PDECHAR), where PDECHAR is a string
        %   outputs two function handles L and N which represent the linear and
        %   the nonlinear parts of the PDE represented by PDECHAR, the time 
        %   interval TSPAN and the initial condition U0.

        % Allen-Cahn equation:
        if ( strcmpi(pdechar, 'AC') == 1 )
            L = @(u) 1e-2*lap(u);
            N = @(u) u - u.^3;
            tspan = [0 60];
            u0 = spherefun(@(x,y,z) cos(cosh(5*x.*z)-10*y));
            
        % Ginzburg-Landau equation:
        elseif ( strcmpi(pdechar, 'GL') == 1 )
            L = @(u) 1e-3*lap(u);
            N = @(u) u - (1 + 1.5i)*u.*(abs(u).^2);
            u0 = randnfunsphere(.1);
            u0 = u0/norm(u0, inf);
            tspan = [0 100];
            
        % Gierer-Meinhardt equations:
        elseif ( strcmpi(pdechar, 'GM') == 1 )
            L = @(u,v) [4e-3*lap(u);4e-2*lap(v)];
            N = @(u,v) [u.^2./v-u;u.^2-v]; 
            tspan = [0 80];
            u01 = @(x,y,z) 1 + .1*(cos(20*x)+cos(20*z)+cos(20*y));
            u02 = @(x,y,z) 1 - .1*(cos(20*x)+cos(20*z)+cos(20*y));
            u0 = chebmatrix([]);
            u0(1,1) = spherefun(u01); 
            u0(2,1) = spherefun(u02);
              
        % Focusing nonlinear Schroedinger equation
        elseif ( strcmpi(pdechar, 'NLS') == 1 )
            L = @(u) 1i*lap(u);
            N = @(u) 1i*u.*abs(u).^2;
            tspan = [0 3];
            A = 1; B = 1;
            u0 = @(lam,th) (2*B^2./(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*th)) - 1)*A;
            u0 = spherefun(u0);
            u0 = .1*u0 + spherefun.sphharm(8, 6);
            u0 = rotate(u0/norm(u0, inf), 0, -pi/4, 0);
        else
            error('SPINOPSPHERE:parseInputs', 'Unrecognized PDE.')
        end

    end