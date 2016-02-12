classdef spinop2 < spinoperator
%SPINOP2   Class for representing the spatial part of 2D differential operators 
%for time-dependent PDEs.
%   SPINOP2 is a class for representing the spatial part S of a time-dependent 
%   PDE of the form u_t = S(u) = Lu + N(u) in 2D, where L is a linear 
%   operator and N is a nonlinear operator. 
%
%   S = SPINOP2(PDECHAR) creates a SPINOP2 object S defined by the string 
%   PDECHAR. The default domain is [-1 1]^2. Strings available include 'GL2' for 
%   Ginzburg-Landau equation and 'GS2' for Gray-Scott equations. Many other PDEs 
%   are available, see HELP/SPIN2.
%
%   S = SPINOP2(PDEFUNLIN, PDEFUNNONLIN) creates a SPINOP2 object S defined by 
%   the function handles PDEFUNLIN (representing the linear part L) and 
%   PDEFUNNONLIN (representing N). The default domain is [-1 1]^2. See Remark 1.
%
%   S = SPINOP2(..., DOM) creates a SPINOP2 object S on the domain DOM.
%
% Remark 1: The linear part PDEFULIN has to be of the form 
%           
%               @(u) A*laplacian(u) + B*biharmonic(u), 
%
%           for some numbers A and B, and the nonlinear part has to be of the 
%           form, 
%
%               @(u) f(u), 
%
%           where f is a nonlinear function of u that does not involve any 
%           derivatives of u.
%
% See also SPINOPERATOR, SPINOP2, SPINOP3, SPIN2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinop2(varargin)
            
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
                S.domain = [-1 1 -1 1];
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
        if ( strcmpi(pdechar, 'GL2') == 1 )
            L = @(u) lap(u);
            N = @(u) u - (1 + 1.3i)*u.*(abs(u).^2);
            dom = [0 200 0 200];
            tspan = [0 80];
            vals = .1*randn(128, 128);
            u0 = chebfun2(vals, dom, 'trig');
            
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS2') == 1 )
            L = @(u,v) [2e-5*lap(u); 1e-5*lap(v)];
            N = @(u,v) [3.5e-2*(1 - u) - u.*v.^2; -9.5e-2*v + u.*v.^2];
            G = 1.25;
            dom = G*[0 1 0 1];
            tspan = [0 3200];
            u01 = @(x,y) 1 - exp(-150*((x-G/2).^2 + (y-G/2).^2));
            u01 = chebfun2(u01, dom, 'trig');
            u02 = @(x,y) exp(-150*((x-G/2).^2 + 2*(y-G/2).^2));
            u02 = chebfun2(u02, dom, 'trig');
            u0 = chebmatrix(u01);
            u0(2,1) = u02;
    
        % Schnakenberg equations:
        elseif ( strcmpi(pdechar, 'Schnak2') == 1 )
            L = @(u,v) [lap(u); 10*lap(v)];
            N = @(u,v) [.1 - u + u.^2.*v; .9 - u.^2.*v];
            G = 50;
            dom = G*[0 1 0 1];
            tspan = [0 200];
            u01 = @(x,y) 1 - exp(-10*((x-G/2).^2 + (y-G/2).^2));
            u01 = chebfun2(u01, dom, 'trig');
            u02 = @(x,y) exp(-10*((x-G/2).^2 + 2*(y-G/2).^2));
            u02 = chebfun2(u02, dom, 'trig');
            u0 = chebmatrix(u01);
            u0(2,1) = u02;
    
        % Swift-Hohenberg equation:
        elseif ( strcmpi(pdechar, 'SH2') == 1 )
            L = @(u) -2*lap(u) - biharm(u);
            N = @(u) -.9*u + u.^2 - u.^3;
            G = 50;
            dom = G*[0 1 0 1];
            tspan = [0 200];
            vals = .1*randn(64, 64);
            u0 = chebfun2(vals, dom, 'trig');
            
        else
            error('SPINOP2:parseInputs', 'Unrecognized PDE.')
        end

    end