classdef spinop < spinoperator
%SPINOP   Class for representing the spatial part of 1D differential operators 
%for time-dependent PDEs.
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
                    S.numVars = nargin(item);
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
            
        % Cahn-Hilliard equation:
        elseif ( strcmpi(pdechar, 'CH') == 1 )
            L = @(u) -0.01*(diff(u, 2) + 0.001*diff(u, 4));
            N = @(u) 0.01*diff(u.^3, 2);
            
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS') == 1 )
            L = @(u,v) [diff(u,2); 0.01*diff(v,2)];
            N = @(u,v) [ 0.02*(1 - u) - u.*v.^2; -.0862*v + u.*v.^2];
            
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
            r = .9; alpha = 10; beta = 2.6;
            L = @(u) (1-r)*diff(u, 2) + alpha*diff(u, 3) + diff(u, 4) + ...
                beta*diff(u, 5) + diff(u, 6);
            N = @(u) -.5*diff(u.^2);
            
       % Nonlinear Schrodinger equation:
        elseif ( strcmpi(pdechar, 'NLS') == 1 )
            L = @(u) 1i*diff(u, 2);
            N = @(u) 1i*abs(u).^2.*u;
            
        % Ohta-Kawasaki equation:
        elseif ( strcmpi(pdechar, 'OK') == 1 )
            ep = .1; sig = 4;
            L = @(u) -ep^2*diff(u, 4) - diff(u, 2) - sig*(u - diff(u, 7));
            N = @(u) diff(u.^3, 2);

        else
            error('SPINOP:getLinearAndNonlinearParts', 'Unrecognized PDE.')
        end

    end