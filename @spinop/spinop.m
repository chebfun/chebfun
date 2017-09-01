classdef spinop < spinoperator
%SPINOP   Class for representing 1D differential operators in SPIN.
%   SPINOP is a class for representing the spatial part S of a time-dependent 
%   PDE of the form u_t = S(u) = Lu + N(u) in 1D, where L is a linear 
%   operator and N is a nonlinear operator. 
%
%   S = SPINOP(PDECHAR) creates a SPINOP object S defined by the string PDECHAR.
%   Strings available include 'AC' for Allen-Cahn equation, 'KS' for 
%   Kuramoto-Sivashinsky equation, and 'KdV' for Korteweg-de Vries equation. 
%   Other PDEs are available, see HELP/SPIN.
%
%   S = SPINOP(DOM, TSPAN) creates a SPINOP object S on DOM x TSPAN. The other
%   fields of a SPINOP are its linear part S.LIN, its nonlienar part S.NONLIN 
%   and the initial condition S.INIT. The fields can be set via S.PROP = VALUE. 
%   See Remark 1 and Example 1.
% 
% Remark 1: The linear part can be any linear constant-coefficient differential
%           operator, e.g., @(u) diff(u,3) + u. The nonlinear part has to be
%           constant-coefficient and of the form @(u) diff(f(u), m) where f is a 
%           nonlinear function of u that does not involve any derivative of u, 
%           and m is the differentiation order. The following function handles 
%           are allowed:
%
%               N = @(u) .5*diff(u.^2);
%               N = @(u) diff(u.^2 + u.^3, 2);
%               N = @(u) exp(u) + cos(sin(2*u));
%
%           The following function handles are not allowed:
%
%               N = @(u) u.*diff(u); 
%               N = @(u) diff(u.^2, 2) + diff(u.^3, 2);
% 
%           For systems of equations, the nonlinear part has to be of the form 
%           @(u) f(u), i.e., no differentiation.         
%          
% Example 1: To construct a SPINOP corresponding to the KdV equation on 
%            DOM = [-pi pi] x TSPAN = [0 1e-2] with initial condition 
%            u0(x) = 100*sin(x), one can type
%
%            dom = [-pi pi]; tspan = [0 1e-2];
%            S = spinop(dom, tspan);
%            S.lin = @(u) -diff(u,3);
%            S.nonlin = @(u) -.5*diff(u.^2);
%            S.init = chebfun(@(x) 100*sin(x), dom);
%
% See also SPINOPERATOR, SPIN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinop(varargin)
            
            % Zero input argument:
            if ( nargin == 0 )
                return
            
            % One input argument:
            elseif ( nargin == 1 )
                if ( isa(varargin{1}, 'char') == 1 )
                    [L, N, dom, tspan, u0] = parseInputs(varargin{1});
                    S.lin = L;
                    S.nonlin = N;
                    S.domain = dom;
                    S.tspan = tspan;
                    S.init = u0;
                else
                    error('SPINOP:constructor', ['When constructing a ', ...
                        'SPINOP with one input argument, this argument ', ...
                        'should be a STRING.'])
                end
            
            % Two input arguments:
            elseif ( nargin == 2 )
                countDouble = 0;
                for j = 1:nargin
                    item =  varargin{j};
                    if ( isa(item, 'double') == 1 && countDouble == 0 )
                        S.domain = item;
                        countDouble = 1;
                    elseif ( isa(item, 'double') == 1 && countDouble == 1 )
                        S.tspan = item;
                    else
                    error('SPINOP:constructor', ['When constructing a ', ...
                        'SPINOP with two input arguments, these arguments ', ...
                        'should be DOUBLE.'])
                    end
                end
            
            % More than two input arguments:
            else
                error('SPINOP:constructor', 'Too many input arguments.')
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
            tspan = [0 500];
            u0 = @(x) 1/3*tanh(2*sin(x)) - exp(-23.5*(x-pi/2).^2) ...
                + exp(-27*(x-4.2).^2) + exp(-38*(x-5.4).^2);
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
            u02 = chebfun(@(x) exp(-100*x.^2), dom, 'trig');
            u03 = chebfun(@(x) exp(-100*(x-.5).^2), dom, 'trig');
            u0 = [u01; u02; u03];
            
        % Cahn-Hilliard equation:
        elseif ( strcmpi(pdechar, 'CH') == 1 )
            L = @(u) -1e-2*(diff(u, 2) + 1e-3*diff(u, 4));
            N = @(u) 1e-2*diff(u.^3, 2);
            dom = [-1 1];
            tspan = [0 100];
            u0 = chebfun('1/5*(sin(4*pi*x)).^5 - 4/5*sin(pi*x)', dom, 'trig');
            
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS') == 1 )
            L = @(u,v) [diff(u,2); 1e-2*diff(v,2)];
            N = @(u,v) [2e-2*(1 - u) - u.*v.^2; -8.62e-2*v + u.*v.^2];
            G = 50;
            dom = G*[-1 1];
            tspan = [0 8000];
            u01 = chebfun(@(x) 1 - 1/2*sin(pi*(x-G)/(2*G)).^100, dom, 'trig');
            u02 = chebfun(@(x) 1/4*sin(pi*(x-G)/(2*G)).^100, dom, 'trig');
            u0 = [u01; u02];
            
        % Korteweg-de Vries equation:
        elseif ( strcmpi(pdechar, 'KdV') == 1 )
            L = @(u) -diff(u, 3);
            N = @(u) -.5*diff(u.^2);
            dom = [-pi pi];
            A = 25; B = 16;
            tspan = [0 0.03015];
            u0 = @(x) 3*A^2*sech(.5*A*(x+2)).^2 + 3*B^2*sech(.5*B*(x+1)).^2;
            u0 = chebfun(u0, dom, 'trig');
            
        % Kuramoto-Sivashinsky equation:
        elseif ( strcmpi(pdechar, 'KS') == 1 )
            L = @(u) -diff(u, 2) - diff(u, 4);
            N = @(u) -.5*diff(u.^2);
            dom = [0 32*pi];
            tspan = [0 200];
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
            tspan = [0 20];
            A = 1; B = 1;
            u0 = @(x) (2*B^2./(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*x)) - 1)*A;
            u0 = chebfun(u0, dom, 'trig');
             
        else
            error('SPINOP:parseInputs', 'Unrecognized PDE.')
        end

    end
