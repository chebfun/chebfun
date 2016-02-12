classdef spinop3 < spinoperator
%SPINOP3   Class for representing the spatial part of 3D differential operators 
%for time-dependent PDEs.
%   SPINOP3 is a class for representing the spatial part S of a time-dependent 
%   PDE of the form u_t = S(u) = Lu + N(u) in 3D, where L is a linear 
%   operator and N is a nonlinear operator. 
%
%   S = SPINOP3(PDECHAR) creates a SPINOP3 object S defined by the string 
%   PDECHAR. The default domain is [-1 1]^3. Strings available include 'GL3' for
%   Ginzburg-Landau equation and 'GS3' for Gray-Scott equations. Many other PDEs
%   are available, see HELP/SPIN3.
%
%   S = SPINOP3(PDEFUNLIN, PDEFUNNONLIN) creates a SPINOP3 object S defined by 
%   the function handles PDEFUNLIN (representing the linear part L) and 
%   PDEFUNNONLIN (representing N). The default domain is [-1 1]^3. See Remark 1.
%
%   S = SPINOP3(..., DOM) creates a SPINOP3 object S on the domain DOM.
%
% Remark 1: The linear part PDEFULIN has to be of the form 
%           
%               @(u) A*laplacian(u) + B*biharmonic(u), 
%
%           for some numbers A and B, and the nonlinear part has to be of the 
%           form, 
%
%           where f is a nonlinear function of u that does not involve any 
%           derivatives of u.
%
% See also SPINOPERATOR, SPINOP, SPINOP2, SPIN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinop3(varargin)
            
            if ( nargin == 0 )
                return
            end
            
            countFun = 0;
            for j = 1:nargin
                item =  varargin{j};
                if ( isa(item, 'char') == 1 )
                    [L, N, dom, tspan, u0] = parseInputs(item);
                    S.linearPart = L;
                    S.nonlinearPart = N;
                    S.domain = dom;
                    S.tspan = tspan;
                    S.init = u0;
                elseif ( isa(item, 'function_handle') && countFun == 0 )
                    S.linearPart = item;
                    countFun = 1;
                elseif ( isa(item, 'function_handle') && countFun == 1 )
                    S.nonlinearPart = item;
                elseif ( isa(item, 'double') == 1 )
                    S.domain = item;
                end
            end
            
            if ( isempty(S.domain) == 1 )
                S.domain = [-1 1 -1 1 -1 1];
            end
            
        end
        
    end
    
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% METHODS IMPLEMENTED IN THIS FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    function [L, N, dom, tspan, u0] = parseInputs(pdechar)
        %PARSEINPUTS   Parse inputs when using a STRING.
        %   [L, N, DOM, TSPAN, U0] = PARSEINPUTS(PDECHAR), where PDECHAR is a 
        %   string, outputs two function handles L and N which represent the 
        %   linear and the nonlinear parts of the PDE represented by PDECHAR, 
        %   the domain DOM, the time interval TSPAN and the initial condition
        %   U0.
            
        % Ginzburg-Landau equation:
        if ( strcmpi(pdechar, 'GL3') == 1 )
            L = @(u) lap(u);
            N = @(u) u - (1 + 1.3i)*u.*(abs(u).^2);
            dom = [0 100 0 100 0 100];
            tspan = [0 50];
            vals = .1*randn(32, 32, 32);
            u0 = chebfun3(vals, dom, 'trig');
            
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS3') == 1 )
            L = @(u,v) [2e-5*lap(u); 1e-5*lap(v)];
            N = @(u,v) [3.5e-2*(1 - u) - u.*v.^2; -9.5e-2*v + u.*v.^2];
            G = 0.75;
            dom = G*[0 1 0 1 0 1];
            tspan = [0 1600];
            u01 = @(x,y,z) 1 - exp(-150*((x-G/2).^2 + (y-G/2).^2 + (z-G/2).^2));
            u01 = chebfun3(u01, dom, 'trig');
            u02 = @(x,y,z) exp(-150*((x-G/2).^2 + 2*(y-G/2).^2 + (z-G/2).^2));
            u02 = chebfun3(u02, dom, 'trig');
            u0 = chebmatrix(u01);
            u0(2,1) = u02;
            
        % Schnakenberg equations:
        elseif ( strcmpi(pdechar, 'Schnak3') == 1 )
            L = @(u,v) [lap(u); 10*lap(v)];
            N = @(u,v) [.1 - u + u.^2.*v; .9 - u.^2.*v];
            G = 20;
            dom = G*[0 1 0 1 0 1];
            tspan = [0 200];
            u01 = @(x,y,z) 1 - exp(-10*((x-G/2).^2 + (y-G/2).^2 + (z-G/2).^2));
            u01 = chebfun3(u01, dom, 'trig');
            u02 = @(x,y,z) exp(-10*((x-G/2).^2 + 2*(y-G/2).^2 + (z-G/2).^2));
            u02 = chebfun3(u02, dom, 'trig');
            u0 = chebmatrix(u01);
            u0(2,1) = u02;
            
        % Swift-Hohenberg equation:
        elseif ( strcmpi(pdechar, 'SH3') == 1 )
            L = @(u) -2*lap(u) - biharm(u);
            N = @(u) -.9*u + u.^2 - u.^3;
            G = 50;
            dom = G*[0 1 0 1 0 1];
            tspan = [0 200];
            u0 = @(x,y,z) 1/4*(sin(pi*x/10) + sin(pi*y/10) + sin(pi*z/10) + ...
                sin(pi*x/2).*sin(pi*y/2) + sin(pi*x/2).*sin(pi*z/2) + ...
                sin(pi*z/2).*sin(pi*y/2));
            u0 = chebfun3(u0, dom);
    
        else
            error('SPINOP3:parseInputs', 'Unrecognized PDE.')
        end

    end