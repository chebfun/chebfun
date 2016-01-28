classdef spinop < spinoperator
%SPINOP   Class for representing the spatial part of 1D differential operators 
%for time-dependent PDEs.
%   SPINOP is a class for representing the spartial part S of a time-dependent 
%   PDE of the form u_t = S(u) = Lu + N(u) in 1D, where L is a linear 
%   operator and N is a nonlinear operator. 
%
%   S = SPINOP(PDECHAR) creates a SPINOP object S defined by the string PDECHAR.
%   Default domain is [-1 1]. Strings available include 'AC' for Allen-Cahn 
%   equation, 'KS' for Kuramoto-Sivashinsky equation, and 'KdV' for Korteweg-de 
%   Vries equation. Many other PDEs are available, see Remark 1.
%
%   S = SPINOP(PDEFUNLIN, PDEFUNNONLIN) creates a SPINOP object S defined by the 
%   function handles PDEFUNLIN (representing the linear part L) and PDEFUNNONLIN
%   (representing N). Default domain is [-1 1]. See Remark 2.
%
%   S = SPINOP(..., DOM) creates a SPINOP object S on the domain DOM.
%
% Remark 1: Available strings PDECHAR are
%
%    - 'AC' for Allen-Cahn equation 
%
%           u_t = 5e-3*u_xx - u^3 + u,
%
%    - 'Burg' for viscous Burgers equation 
%
%           u_t = 1e-3*u_xx - u*u_x,
%
%    - 'BZ' for Belousov-Zhabotinsky equation
%
%           u_t = 1e-5*diff(u,2) + u + v - u*v - u^2,
%           v_t = 2e-5*diff(v,2) + w - v - u*v,
%           w_t = 1e-5*diff(w,2) + u - w
% 
%    - 'CH' for Cahn-Hilliard equation 
%
%           u_t = 1e-2*(-u_xx - 1e-3*u_xxxx + (u^3)_xx),
%
%    - 'GS' for Gray-Scott equations
%
%           u_t = diff(u,2) + 2e-2*(1-u) - u*v^2,
%           v_t = 1e-2*diff(u,2) - 8.62e-2*v + u*v^2,
%           
%    - 'KdV' for Korteweg-de Vries equation 
%
%           u_t = -u*u_x - u_xxx,
%
%    - 'KS' for Kuramoto-Sivashinsky equation 
%
%           u_t = -u*u_x - u_xx - u_xxxx,
% 
%    - 'Niko' for Nikolaevskiy equation 
%
%           u_t = -u*u_x + 1e-1*u_xx + u_xxxx + u_xxxxxx,
%
%    - 'NLS' for the focusing nonlinear Schrodinger equation 
%
%           u_t = 1i*u_xx + 1i*|u|^2*u,
%
%    - and 'OK' for the Ohta-Kawasaki equation
%
%           u_t = -u_xx - 1e-2*u_xxxx - 4(u - sum(u)) + (u^3)_xx.
%
% Remark 2: The linear part PDEFUNLIN can be any linear constant-coefficient 
%           differential operator, e.g., @(u) diff(u,3) + 5*u or @(u) diff(u,6).
%           The nonlinear part PDEFUNNONLIN has to be constant-coefficient and
%           of the form @(u) C*diff(f(u), m), where C is a number, f is a 
%           nonlinear function of u that does not involve any derivative of u, 
%           and m is the differentiation order. If m=0, it can be of any form. 
%           The following function handles are allowed:
%
%               pdefunnonlin = @(u) .5*diff(u.^2);
%               pdefunnonlin = @(u) diff(u.^2 + u.^3, 2);
%               pdefunnonlin = @(u) exp(u) + cos(sin(2*u));
%
%           The following function handles are not:
%
%               pdefunnonlin = @(u) u.*diff(u); 
%               pdefunnonlin = @(u) diff(u.^2, 2) + diff(u.^3, 2);
% 
%           For example, to construct a SPINOP which correponds to the KdV 
%           equation on [-1 1],       
%           
%               L = @(u) -diff(u, 3);
%               N = @(u) -.5*diff(u.^2);
%               S = spinop(L, N);
% 
%           For systems of equations, PDEFUNNONLIN has to be of the form 
%           @(u) f(u), i.e., no differentiation.         
%          
% See also SPINOPERATOR, SPINOP2, SPINOP3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinop(varargin)
            
            S.dimension = 1;
            
            if ( nargin == 0 )
                return
            end
            
            countFun = 0;
            for j = 1:nargin
                item =  varargin{j};
                if ( isa(item, 'char') == 1 )
                    [L, N] = getLinearAndNonlinearParts(item);
                    S.linearPart = L;
                    S.nonlinearPart = N;
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

    function [L, N] = getLinearAndNonlinearParts(pdechar)
        %GETLINEARANDNONLINEARPARTS   Get the linear and the nonlinear parts of 
        %a 1D time-dependent PDE from a key word.
        %   [L, N] = GETLINEARANDNONLINEARPARTS(PDECHAR), where PDECHAR is a 
        %   string, outputs two function handles L and N which represent the
        %   linear and the nonlinear parts of the PDE represented by PDECHAR.

        % Allen-Cahn equation:
        if ( strcmpi(pdechar, 'AC') == 1 )
            L = @(u) 5e-3*diff(u, 2);
            N = @(u) u - u.^3;
            
        % Viscous Burgers equation:
        elseif ( strcmpi(pdechar, 'Burg') == 1 )
            L = @(u) 1e-3*diff(u, 2);
            N = @(u) -.5*diff(u.^2);
            
        % Belousov-Zhabotinsky equation:
        elseif ( strcmpi(pdechar, 'BZ') == 1 )
            L = @(u,v,w)[1e-5*diff(u,2); 2e-5*diff(v,2); 1e-5*diff(w,2)];
            N = @(u,v,w)[u + v - u.*v - u.^2; w - v - u.*v; u - w];
            
        % Cahn-Hilliard equation:
        elseif ( strcmpi(pdechar, 'CH') == 1 )
            L = @(u) -1e-2*(diff(u, 2) + 1e-3*diff(u, 4));
            N = @(u) 1e-2*diff(u.^3, 2);
            
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS') == 1 )
            L = @(u,v) [diff(u,2); 1e-2*diff(v,2)];
            N = @(u,v) [2e-2*(1 - u) - u.*v.^2; -8.62e-2*v + u.*v.^2];
            
        % Korteweg-de Vries equation:
        elseif ( strcmpi(pdechar, 'KdV') == 1 )
            L = @(u) -diff(u, 3);
            N = @(u) -.5*diff(u.^2);
            
        % Kuramoto-Sivashinsky equation:
        elseif ( strcmpi(pdechar, 'KS') == 1 )
            L = @(u) -diff(u, 2) - diff(u, 4);
            N = @(u) -.5*diff(u.^2);
            
        % Nikolaevskiy equation:
        elseif ( strcmpi(pdechar, 'Niko') == 1 )
            L = @(u) .1*diff(u, 2) + diff(u, 4) + diff(u, 6);
            N = @(u) -.5*diff(u.^2);
            
       % Nonlinear Schrodinger equation:
        elseif ( strcmpi(pdechar, 'NLS') == 1 )
            L = @(u) 1i*diff(u, 2);
            N = @(u) 1i*abs(u).^2.*u;
            
        % Ohta-Kawasaki equation:
        elseif ( strcmpi(pdechar, 'OK') == 1 )
            L = @(u) - diff(u, 2) - 1e-2*diff(u, 4)  - 4*(u - diff(u, 7));
            N = @(u) diff(u.^3, 2);

        else
            error('SPINOP:getLinearAndNonlinearParts', 'Unrecognized PDE.')
        end

    end