classdef spinoperator
%SPINOPERATOR   Abstract class for representing the spatial part of differential 
%operators for time-dependent PDEs. 
%   SPINOPERATOR is a class for representing the spartial part S of a 
%   time-dependent PDE of the form u_t = S(u) = Lu + N(u), where L is a linear 
%   operator and N is a nonlinear operator. SPINOP (in 1D), SPINOP2 (in 2D),
%   SPINOP3 (in 3D) and SPINOPSPHERE (on the sphere) are full implementations.
%   
% See also SPINOP, SPINOP2, SPINOP3, SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        domain              % Spatial domain (1x2 DOUBLE in 1D, 1x4 in 2D and on
                            % the sphere, 1x6 in 3D)
        init                % Initial condition (CHEBFUN in 1D, CHEBFUN2 in 2D,
                            % CHEBFUN3 in 3D, SPHEREFUN on the SPHERE, 
                            % CHEBMATRIX for systems)
        lin                 % Linear part of the operator (FUNCTION HANDLE)
        nonlin              % Nonlinear part of the operator (FUNCTION HANDLE)
        tspan               % Vector of strictly increasing time 
                            % samples (DOUBLE)
    end
    
    % DEPENDENT PROPERTIES:
    properties ( Access = public, Dependent = true )
        numVars             % Number of unknown functions (1x1 INT)
    end
    
    % DEPENDENT AND HIDDEN PROPERTIES:
    properties ( Access = public, Hidden = true, Dependent = true )
        nonlinearPartCoeffs % Differential part of the operator in coefficient
                            % space (FUNCTION HANDLE)
        nonlinearPartVals   % Nondifferential part of the operator in value
                            % space (FUNCTION HANDLE)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS FOR DEPENDENT PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % METHOD for numVars:
        function nVars = get.numVars(S)
            
            nVars = nargin(S.lin);
            
        end
        
        % METHOD for nonlinearPartCoeffs:
        function Nc = get.nonlinearPartCoeffs(S)
            
            % Extract the nonlinear part:
            N = S.nonlin;
            
            % If N is empty, we are done:
            if ( isempty(N) == 1 )
                Nc = @(u) 0*u;
                return
            end
            
            % Else, get the variables of the workspace:
            func = functions(N);
            wrk = func.workspace{1};
            names = fieldnames(wrk);
            if ( isempty(names) == 0 )
                lengthNames = size(names, 1);
                for iter = 1:lengthNames
                    eval(sprintf('%s = wrk.(names{iter});', names{iter}));
                end
            end
            
            % Convert FUNCN to a string STRN:
            strN = func2str(N);
            
            % Get the number of variables NVARS:
            nVars = nargin(N);
            
            % For scalar equations in 1D, we support nonlinearities of the form
            % diff(f(u),m) with m>=0:
            if ( isa(S, 'spinop') == 1 && nVars == 1 )
                
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
                
            % For scalar equations in 2D/3D and on the sphere, and for systems 
            % of equations in 1D/2D/3D and on the sphere, we only support 
            % nonlinearities of the form f_i(u_1,...,u_n), i.e., with no 
            % differentiation, so Nc=1:
            else
                
                % Nc=1:
                Nc = 1;
                
            end
            
        end
        
        % METHOD for nonlinearPartVals:
        function Nv = get.nonlinearPartVals(S)
            
            % Extract the nonlinear part:
            N = S.nonlin;
            
            % If N is empty, we are done:
            if ( isempty(N) == 1 )
                Nv = @(u) 0*u;
                return
            end
            
            % Else, get the variables of the workspace:
            func = functions(N);
            wrk = func.workspace{1};
            names = fieldnames(wrk);
            if ( isempty(names) == 0 )
                lengthNames = size(names, 1);
                for iter = 1:lengthNames
                    eval(sprintf('%s = wrk.(names{iter});', names{iter}));
                end
            end
            
            % Convert FUNCN to a string STRN:
            strN = func2str(N);
            
            % Get the number of variables NVARS:
            nVars = nargin(N);
            
            % For scalar equations in 1D, we support nonlinearities of the form
            % diff(f(u),m) with m>=0:
            if ( isa(S, 'spinop') == 1 && nVars == 1 )
                
                % Get rid of the differentiation part in STRN to get NV=f(u):
                oldString = {'diff', ',\d*)'};
                newString = {'', ')'};
                Nv = regexprep(strN, oldString, newString);
                Nv = eval(Nv);
                
            % For scalar equations in 2D/3D and on the sphere, and for systems 
            % of equations in 1D/2D/3D and on the sphere, we only support 
            % nonlinearities of the form f_i(u_1,...,u_n), i.e., with no 
            % differentiation, so Nv=N:
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
                for iter = 1:nVars
                    idx1 = [num2str(rats((iter-1)/nVars)), '*', 'length(u)', '+1'];
                    idx2 = [num2str(rats(iter/nVars)), '*', 'length(u)'];
                    strNvNew = strrep(strNv, variablesNames{iter}, ...
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
        
        % Discretize a SPINOPERATOR:
        [L, Nc] = discretize(S, N)

        % Returns indexes for dealiasing procedure (2/3-rule):
        idx = getDealiasingIndexes(S, N, nVars)
        
        % Returns the type of CHEBFUN on which the SPINOPERATOR acts:
        out = getChebfunType(S)
        
        % Returns the transform coeffs -> values:
        F = getCoeffs2ValsTransform(S)
        
        % Returns the spatial dimension:
        dim = getDimension(S)
        
        % Returns a grid correspoding to a SPINOPRERATOR object:
        grid = getGrid(S, N, dom)
        
        % Returns the adequate SPINPREFERENCE object:
        pref = getPreference(S)
        
        % Returns the transform values -> coeffs:
        F = getVals2CoeffsTransform(S)

        % Returns 1 if the linear part of the SPINOPERATOR is diagonal, 
        % 0 otherwise:
        out = isDiag(S)
        
        % Initialize a movie when solving a PDE specified by a SPINOPERATOR:
        [p, opts] = initializeMovie(S, dt, pref, v, compGrid, plotGrid)
        
        % Update the movie when solving a PDE specified by a SPINOPERATOR:
        opts = updateMovie(S, dt, p, options, t, v, compGrid, plotGrid)
         
        % Reshape the data that will be used for constructing the solution at
        % the end of the time-stepping: (For example, for SPINOPSPEHRE, it 
        % extracts half of the data, since the data has been doubled-up with
        % the DFS method.)
        data = reshapeData(S, data, nVars)
        
        % Add the (repeated) endpoints to a periodic grid:
        grid = reshapeGrid(S, grid)

    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE AND STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = true )
        
        % Solve a PDE defined by a SPINOPERATOR:
        [uout, tout, computingTime] = solvepde(varargin)
   
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
