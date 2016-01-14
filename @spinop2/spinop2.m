classdef spinop2 < spinoperator
%SPINOP2   Class for representing the spatial part of 2D differential operators 
%for time-dependent PDEs.
%
% See also SPINOPERATOR, SPINOP, SPINOP3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinop2(varargin)
            
            S.dimension = 2;
            
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
                    S.numVars = nargin(L);
                elseif ( isa(item, 'function_handle') && countFun == 0 )
                    S.linearPart = item;
                    S.numVars = nargin(item);
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

    function [L, N] = getLinearAndNonlinearParts(pdechar)
        %GETLINEARANDNONLINEARPARTS   Get the linear and the nonlinear parts of 
        %a 2D time-dependent PDE from a key word.
        %   [L, N] = GETLINEARANDNONLINEARPARTS(PDECHAR), where PDECHAR is a 
        %   string, outputs two function handles L and N which represent the
        %   linear and the nonlinear parts of the PDE represented by PDECHAR.
            
        % Ginzburg-Landau equations:
        if ( strcmpi(pdechar, 'GL2') == 1 )
            L = @(u) laplacian(u);
            N = @(u) u - (1 + 1.3i)*u.*(abs(u).^2);
        
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS2') == 1 )
            L = @(u,v) [2e-5*laplacian(u); 1e-5*laplacian(v)];
            % Mitosis (1):     F = 0.0367; K = 0.0649;
            % Fingerpint (1):  F = 0.0545; K = 0.062;
            % Mitosis (2):     F = 0.035;  K = 0.0625;
            % Fingerpint (2):  F = 0.035;  K = 0.06;
            F = 0.035;  K = 0.06;
            N = @(u,v) [ F*(1 - u) - u.*v.^2; -(F+K)*v + u.*v.^2];
            
        % Schnakenberg equations:
        elseif ( strcmpi(pdechar, 'Schnak2') == 1 )
            L = @(u,v) [laplacian(u); 10*laplacian(v)];
            A = .1; B = .9; G = 1;
            N = @(u,v) [ G*(A - u + u.^2.*v); G*(B - u.^2.*v)];
            
        else
            error('SPINOP2:getLinearAndNonlinearParts', 'Unrecognized PDE.')
        end

    end