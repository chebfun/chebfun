classdef spinop3 < spinoperator
%SPINOP3   Class for representing the spatial part of 3D differential operators 
%for time-dependent PDEs.
%   SPINOP3 is a class for representing the spartial part S of a time-dependent 
%   PDE of the form u_t = S(u) = Lu + N(u) in 3D, where L is a linear 
%   operator and N is a nonlinear operator. 
%
%   S = SPINOP3(PDECHAR) creates a SPINOP3 object S defined by the string 
%   PDECHAR. Default domain is [-1 1 -1 1 -1 1]. Strings available include 'GL3' 
%   for Ginzburg-Landau equation and 'GS3' for Gray-Scott equations. Many other 
%   PDEs are available, see Remark 1.
%
%   S = SPINOP3(PDEFUNLIN, PDEFUNNONLIN) creates a SPINOP3 object S defined by 
%   the function handles PDEFUNLIN (representing the linear part L) and 
%   PDEFUNNONLIN (representing N). Default domain is [-1 1 -1 1 -1 1]. See 
%   Remark 2.
%
%   S = SPINOP3(..., DOM) creates a SPINOP3 object S on the domain DOM.
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
% Remark 2: The linear part PDEFULIN has to be of the form 
%           
%               @(u) A*laplacian(u) + B*biharmonic(u), 
%
%           for some numbers A and B, and the linear part has be of the form, 
%
%               @(u) f(u), 
%
%           where f is a nonlinear unction of u that does not involve any 
%           derivative of u.
%
% See also SPINOPERATOR, SPINOP, SPINOP2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinop3(varargin)
            
            S.dimension = 3;
            
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
                S.domain = [-1 1 -1 1 -1 1];
            end
            
        end
        
    end
    
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% METHODS IMPLEMENTED IN THIS FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    function [L, N] = getLinearAndNonlinearParts(pdechar)
        %GETLINEARANDNONLINEARPARTS   Get the linear and the nonlinear parts of 
        %a 3D time-dependent PDE from a key word.
        %   [L, N] = GETLINEARANDNONLINEARPARTS(PDECHAR), where PDECHAR is a 
        %   string, outputs two function handles L and N which represent the
        %   linear and the nonlinear parts of the PDE represented by PDECHAR.
            
        % Ginzburg-Landau equation:
        if ( strcmpi(pdechar, 'GL3') == 1 )
            L = @(u) laplacian(u);
            N = @(u) u - (1 + 1.3i)*u.*(abs(u).^2);
        
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS3') == 1 )
            L = @(u,v) [2e-5*laplacian(u); 1e-5*laplacian(v)];
            N = @(u,v) [3.5e-2*(1 - u) - u.*v.^2; -9.5e-2*v + u.*v.^2];
            
        % Schnakenberg equations:
        elseif ( strcmpi(pdechar, 'Schnak3') == 1 )
            L = @(u,v) [laplacian(u); 10*laplacian(v)];
            N = @(u,v) [.1 - u + u.^2.*v; .9 - u.^2.*v];
                          
        % Swift-Hohenberg equation:
        elseif ( strcmpi(pdechar, 'SH3') == 1 )
            L = @(u) -2*laplacian(u) - biharmonic(u);
            N = @(u) -.9*u + u.^2 - u.^3;
            
        else
            error('SPINOP3:getLinearAndNonlinearParts', 'Unrecognized PDE.')
        end

    end