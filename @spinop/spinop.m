classdef spinop < spinoperator
%SPINOP   Class for representing the spatial part of 1D differential operators 
%for time-dependent PDEs.
%
% See also SPINOPERATOR, SPINOP2, SPINOP3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        dimension = 1; % Spatial dimension (1 for SPINOP)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinop(varargin)

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
            
            [Nc, Nv] = getNonlinearPartsCoeffsAndVals(S.nonlinearPart);
            S.nonlinearPartCoeffs = Nc;
            S.nonlinearPartVals = Nv;
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE AND NON-STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = false )
        
    end

end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% METHODS IMPLEMENTED IN THIS FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

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
            
        % Cahn-Hilliard equation:
        elseif ( strcmpi(pdechar, 'CH') == 1 )
            L = @(u) -0.01*(diff(u, 2) + 0.001*diff(u, 4));
            N = @(u) 0.01*diff(u.^3, 2);
            
        % Gray-Scott equations:
        elseif ( strcmpi(pdechar, 'GS') == 1 )
            L = @(u,v) [diff(u,2); 0.01*diff(v,2)];
            N = @(u,v) [ 0.09*(1 - u) - u.*v.^2; - 0.0862*v + u.*v.^2];
            
        % Korteweg-de Vries equation:
        elseif ( strcmpi(pdechar, 'KdV') == 1 )
            L = @(u) -diff(u, 3);
            N = @(u) -.5*diff(u.^2);
            
        % Kuramoto-Sivashinsky equation:
        elseif ( strcmpi(pdechar, 'KS') == 1 )
            L = @(u) -(diff(u, 2) + diff(u, 4));
            N = @(u) -.5*diff(u.^2);
            
       % Nonlinear Schrodinger equation:
        elseif ( strcmpi(pdechar, 'NLS') == 1 )
            L = @(u) 1i*diff(u, 2);
            N = @(u) 1i*abs(u).^2.*u;
 
        % Viscous Burgers equation:
        elseif ( strcmpi(pdechar, 'Burg') == 1 )
            L = @(u) 1e-3*diff(u, 2);
            N = @(u) -.5*diff(u.^2);

        else
            error('SPINOP:getLinearPart', 'Unrecognized PDE.')
        end

    end
    
    function [Nc, Nv] = getNonlinearPartsCoeffsAndVals(N)
        %GETLINEARANDNONLINEARPARTS   Get the nonlinear parts in coefficient and
        %value space of a 1D time-dependent PDE.
        %   [Nc, Nv] = GETNONLINEARPARTSCOEFFSANDVALS(N), where N is a function
        %   handle, outputs two function handles NC and NV which represent 
        %   the nonlinear parts, in coefficient and value space, of the 
        %   nonlienar part of the PDE represented by N.
        
        % Get the variables of the workspace:
        func = functions(N);
        wrk = func.workspace{1};
        names = fieldnames(wrk);
        if ( isempty(names) == 0 )
            lengthNames = size(names, 1);
            for k = 1:lengthNames
                eval(sprintf('%s = wrk.(names{k});', names{k}));
            end
        end
        
        % Convert FUNCN to a string STRN:
        strN = func2str(N);
        
        % Get the number of variables NVARS:
        nVars = nargin(N);
        
        % For scalar equation, we support nonlinearities of the form
        %  diff(f(u),m) with m>=0:
        if ( nVars == 1)
            
            % Get rid of the differentiation part in STRN to get NV=f(u):
            oldString = {'diff', ',\s*\d*)'};
            newString = {'', ')'};
            Nv = regexprep(strN, oldString, newString);
            Nv = str2func(Nv);
            
            % Compute the differentiation order to get NC=diff(u,m):
            C = chebop(str2func(strrep(strN, '@(', '@(x,')));
            L = linearize(C);
            m = getDiffOrder(L);
            Nc = ['@(u) diff(u,', num2str(m), ')'];
            Nc = str2func(Nc);
        
        % For systems of equations, we only support nonlinearities of the form
        % f_i(u_1,...,u_n), i.e., with no differentiation, so Nv=N and Nc=1:
        else
            
            % Nv=N and Nc=1:
            strNv = func2str(N);      
            Nc = 1;
            
            % We're going to relabel the variables, e.g., @(u,v) u + v is going
            % to be relabelled @(u) u(1:length(u)/2) + u(length(u)/2+1:end).
            % First, get the names of the variables:
            openParenthesis = strfind(strNv, '(');
            openParenthesis = openParenthesis(1);
            closeParenthesis = strfind(strNv, ')');
            closeParenthesis = closeParenthesis(1);
            variablesNames = strNv(openParenthesis+1:closeParenthesis-1);
            variablesNames = regexp(variablesNames,  ',', 'split');
            
            % Second, relabel the variables:
            strNv = strNv(closeParenthesis+1:end);
            for k = 1:nVars
                idx1 = [num2str((k-1)/nVars), '*', 'length(u)', '+1'];
                idx2 = [num2str(k/nVars), '*', 'length(u)'];
                strNvNew = strrep(strNv, variablesNames{k}, ...
                    ['u(', idx1, ':', idx2, ')']);
                strNv = strNvNew;
            end
            Nv = str2func(['@(u)', strNvNew]);

        end

    end