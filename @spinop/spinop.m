classdef spinop < spinoperator
%SPINOP   Class for representing the spatial part of 1D differential operators 
%for time-dependent PDEs.
%   SPINOP is a class for representing the spatial part S of a time-dependent 
%   PDE of the form u_t = S(u) = Lu + N(u) in 1D, where L is a linear 
%   operator and N is a nonlinear operator. 
%
%   S = SPINOP(PDECHAR) creates a SPINOP object S defined by the string PDECHAR.
%   The default domain is [-1 1]. Strings available include 'AC' for Allen-Cahn 
%   equation, 'KS' for Kuramoto-Sivashinsky equation, and 'KdV' for Korteweg-de 
%   Vries equation. Many other PDEs are available, see HELP/SPIN.
%
%   S = SPINOP(PDEFUNLIN, PDEFUNNONLIN) creates a SPINOP object S defined by the 
%   function handles PDEFUNLIN (representing the linear part L) and PDEFUNNONLIN
%   (representing N). The default domain is [-1 1]. See Remark 1.
%
%   S = SPINOP(..., DOM) creates a SPINOP object S on the domain DOM.
%
% Remark 1: The linear part PDEFUNLIN can be any linear constant-coefficient 
%           differential operator, e.g., @(u) diff(u,3) + 5*u or @(u) diff(u,6).
%           The nonlinear part PDEFUNNONLIN has to be constant-coefficient and
%           of the form @(u) C*diff(f(u), m), where C is a number, f is a 
%           nonlinear function of u that does not involve any derivatives of u, 
%           and m is the differentiation order. If m=0, PDEFUNNONLIN can be of 
%           any form. The following function handles are allowed:
%
%               pdefunnonlin = @(u) .5*diff(u.^2);
%               pdefunnonlin = @(u) diff(u.^2 + u.^3, 2);
%               pdefunnonlin = @(u) exp(u) + cos(sin(2*u));
%
%           The following function handles are not allowed:
%
%               pdefunnonlin = @(u) u.*diff(u); 
%               pdefunnonlin = @(u) diff(u.^2, 2) + diff(u.^3, 2);
% 
%           For systems of equations, PDEFUNNONLIN has to be of the form 
%           @(u) f(u), i.e., no differentiation.         
%          
% See also SPINOPERATOR, SPINOP2, SPINOP3, SPIN.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinop(varargin)
            
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
                S.domain = [-1 1];
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

        % Allen-Cahn equation:
        if ( strcmpi(pdechar, 'AC') == 1 )
            L = @(u) 5e-3*diff(u, 2);
            N = @(u) u - u.^3;
            dom = [0 2*pi];
            tspan = [0 300];
            u0 = @(x) tanh(2*sin(x)) + 3*exp(-27.*(x-4.2).^2) ...
                - 3*exp(-23.5.*(x-pi/2).^2) + 3*exp(-38.*(x-5.4).^2);
            u0 = chebfun(u0, dom, 'trig');
            
        % Viscous Burgers equation:
        elseif ( strcmpi(pdechar, 'Burg') == 1 )
            L = @(u) 1e-3*diff(u, 2);
            N = @(u) -.5*diff(u.^2);
            dom = [-1 1];
            tspan = [0 20];
            u0 = chebfun('(1-x.^2).*exp(-30.*(x+1/2).^2)', dom, 'trig');
            
        % Belousov-Zhabotinsky equation:
        elseif ( strcmpi(pdechar, 'BZ') == 1 )
            L = @(u,v,w)[1e-5*diff(u,2); 2e-5*diff(v,2); 1e-5*diff(w,2)];
            N = @(u,v,w)[u + v - u.*v - u.^2; w - v - u.*v; u - w];
            dom = [-1 1];
            tspan = [0 30];
            u01 = chebfun(@(x) exp(-100*(x+.5).^2), dom, 'trig');
            u02 = chebfun(@(x) exp(-100*(x).^2), dom, 'trig');
            u03 = chebfun(@(x) exp(-100*(x-.5).^2), dom, 'trig');
            u0 = [u01; u02; u03];
            
        % Cahn-Hilliard equation:
        elseif ( strcmpi(pdechar, 'CH') == 1 )
            L = @(u) -1e-2*(diff(u, 2) + 1e-3*diff(u, 4));
            N = @(u) 1e-2*diff(u.^3, 2);
            dom = [-1 1];
            tspan = [0 70];
            u0 = chebfun('(sin(4*pi*x)).^5 - sin(pi*x)', dom, 'trig');
            
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS') == 1 )
            L = @(u,v) [diff(u,2); 1e-2*diff(v,2)];
            N = @(u,v) [2e-2*(1 - u) - u.*v.^2; -8.62e-2*v + u.*v.^2];
            G = 50;
            dom = G*[-1 1];
            tspan = [0 12000];
            u01 = chebfun(@(x) 1 - 1/2*sin(pi*(x-G)/(2*G)).^100, dom, 'trig');
            u02 = chebfun(@(x) 1/4*sin(pi*(x-G)/(2*G)).^100, dom, 'trig');
            u0 = [u01; u02];
            
        % Korteweg-de Vries equation:
        elseif ( strcmpi(pdechar, 'KdV') == 1 )
            L = @(u) -diff(u, 3);
            N = @(u) -.5*diff(u.^2);
            dom = [-pi pi];
            A = 25^2; B = 16^2;
            tspan = [0 0.03015];
            u0 = @(x) 3*A*sech(.5*sqrt(A)*x).^2 + 3*B*sech(.5*sqrt(B)*(x-1)).^2;
            u0 = chebfun(u0, dom, 'trig');
            
        % Kuramoto-Sivashinsky equation:
        elseif ( strcmpi(pdechar, 'KS') == 1 )
            L = @(u) -diff(u, 2) - diff(u, 4);
            N = @(u) -.5*diff(u.^2);
            dom = [0 32*pi];
            tspan = [0 300];
            u0 = chebfun('cos(x/16).*(1 + sin((x-1)/16))', dom, 'trig');
            
        % Nikolaevskiy equation:
        elseif ( strcmpi(pdechar, 'Niko') == 1 )
            L = @(u) .1*diff(u, 2) + diff(u, 4) + diff(u, 6);
            N = @(u) -.5*diff(u.^2);
            dom = [0 32*pi];
            tspan = [0 300];
            u0 = chebfun('cos(x/16).*(1 + sin((x-1)/16))', dom, 'trig');
            
       % Nonlinear Schroedinger equation:
        elseif ( strcmpi(pdechar, 'NLS') == 1 )
            L = @(u) 1i*diff(u, 2);
            N = @(u) 1i*abs(u).^2.*u;
            dom = [-pi pi];
            tspan = [0 18];
            A = 2; B = 1;
            u0 = @(x) (2*B^2./(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*x)) - 1)*A;
            u0 = chebfun(u0, dom, 'trig');
            
        % Ohta-Kawasaki equation:
        elseif ( strcmpi(pdechar, 'OK') == 1 )
            L = @(u) - diff(u, 2) - 1e-2*diff(u, 4)  - 4*(u - diff(u, 7));
            N = @(u) diff(u.^3, 2);
            dom = [0 2*pi];
            tspan = [0 4];
            u0 = chebfun('cos(x)/2', dom, 'trig');
            
        else
            error('SPINOP:parseInputs', 'Unrecognized PDE.')
        end

    end