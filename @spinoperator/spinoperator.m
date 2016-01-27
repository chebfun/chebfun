classdef spinoperator
%SPINOPERATOR   Abstract class for representing the spatial part of differential 
%operators for time-dependent PDEs. 
%   SPINOPERATOR is a class for representing the spartial part S of a 
%   time-dependent PDE of the form u_t = S(u) = Lu + N(u), where L is a linear 
%   operator and N is a nonlinear operator. SPINOP (in 1D), SPINOP2 (in 2D) and 
%   SPINOP3 (in 3D) are full implementations.
%   
% See also SPINOP, SPINOP2, SPINOP3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        dimension           % Spatial dimension
        domain              % Spatial domain
        linearPart          % Linear part of the operator
        nonlinearPart       % Nonlinear part of the operator
    end
    
    % DEPENDENT PROPERTIES:
    properties ( Access = public, Dependent = true )
        numVars             % Number of unknown functions (>1 for systems)
    end
    
    % DEPENDENT AND HIDDEN PROPERTIES:
    properties ( Access = public, Hidden = true, Dependent = true )
        nonlinearPartCoeffs % Nonlinear part of the operator in coeff. space
        nonlinearPartVals   % Nonlinear part of the operator in value space
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS FOR DEPENDENT PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % METHOD for numVars:
        function nVars = get.numVars(S)
            
            nVars = nargin(S.linearPart);
            
        end
        
        % METHOD for nonlinearPartCoeffs:
        function Nc = get.nonlinearPartCoeffs(S)
            
            % Get the variables of the workspace:
            N = S.nonlinearPart;
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
            
            % Get the dimension DIM:
            dim = S.dimension;
            
            % For scalar equations in 1D, we support nonlinearities of the form
            % diff(f(u),m) with m>=0:
            if ( dim == 1 && nVars == 1 )
                
                % Compute the differentiation order to get NC=diff(u,m):
                diffOrderTwoOrGreater = regexp(strN, ',\d*', 'match');
                diffOrderOne = regexp(strN, 'diff(', 'match');
                if ( isempty(diffOrderTwoOrGreater) == 0 )
                    diffOrder = diffOrderTwoOrGreater{1}(2:end);
                elseif ( isempty(diffOrderOne) == 0 )
                    diffOrder = '1';
                else
                    diffOrder = '0';
                end
                Nc = ['@(u) diff(u,', diffOrder, ')'];
                Nc = eval(Nc);
                
            % For scalar equations in 2D and 3D, and for systems of equations in 
            % 1D, 2D and 3D, we only support nonlinearities of the form 
            % f_i(u_1,...,u_n), i.e., with no differentiation, so Nc=1:
            else
                
                % Nc=1:
                Nc = 1;
                
            end
            
        end
        
        % METHOD for nonlinearPartVals:
        function Nv = get.nonlinearPartVals(S)
            
            % Get the variables of the workspace:
            N = S.nonlinearPart;
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
            
            % Get the dimension DIM:
            dim = S.dimension;
            
            % For scalar equations in 1D, we support nonlinearities of the form
            % diff(f(u),m) with m>=0:
            if ( dim == 1 && nVars == 1 )
                
                % Get rid of the differentiation part in STRN to get NV=f(u):
                oldString = {'diff', ',\d*)'};
                newString = {'', ')'};
                Nv = regexprep(strN, oldString, newString);
                Nv = eval(Nv);
                
            % For scalar equations in 2D and 3D, and for systems of equations in 
            % 1D, 2D and 3D, we only support nonlinearities of the form 
            % f_i(u_1,...,u_n), i.e., with no differentiation, so Nv=N:
            else
                
                % Nv=N:
                strNv = func2str(N);
                
                % We're going to relabel the variables, e.g., @(u,v) u + v is 
                % going to be relabelled @(u) u(1:length(u)/2) + ...
                % u(length(u)/2+1:end). First, get the names of the variables:
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
                        ['u(', idx1, ':', idx2, ',:,:)']);
                    strNv = strNvNew;
                end
                Nv = ['@(u)', strNvNew];
                Nv = eval(Nv);
                
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT AND NON-STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = true, Static = false )
        
        % Check if the solution is resolved in space:
        ishappy = checkHappiness(S, c, pref)
        
        % Discretize a SPINOPERATOR:
        [L, Nc] = discretize(S, N)
        
        % Initialize a movie when solving a PDE specified by a SPINOPERATOR:
        [p, plotOptions] = initializeMovie(S, dt, pref, v, gridPoints)
        
        % Plot a movie when solving a PDE specified by a SPINOPERATOR:
        plotOptions = plotMovie(S, dt, p, plotOptions, t, v, gridPoints)
        
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE AND STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = true )
        
        % Solve a PDE defined by a SPINOPERATOR:
        [uout, tout] = solvepde(varargin)
   
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE AND NON-STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = false )
        
        % Create a contour around each eigenvalue of the linear part of a 
        % SPINOPERATOR:
        LR = computeLR(S, dt, L, M)
           
    end
    
end