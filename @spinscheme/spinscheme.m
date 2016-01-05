classdef spinscheme
%SPINSCHEME   Class for representing timestepping schemes.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        
        internalStages  % number of internal stages
        order           % order of the method
        scheme          % Timestepping scheme
        steps           % number of previous timesteps used (1 if one-step
                        % method, > 1 if multistep method)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function K = spinscheme(schemeName)

            if ( nargin == 0 )
                return
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ETD MULTISTEP:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ETD RUNGE-KUTTA:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ( strcmpi(schemeName, 'etdrk4') == 1 )
                K.order = 4;
                K.internalStages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'exprk5s8') == 1 )
                K.order = 5;
                K.internalStages = 8;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'friedli') == 1 )
                K.order = 4;
                K.internalStages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'hochbruck-ostermann') == 1 )
                K.order = 4;
                K.internalStages = 5;
                K.steps = 1;
                K.scheme = schemeName;

            elseif ( strcmpi(schemeName, 'krogstad') == 1 )
                K.order = 4;
                K.internalStages = 4;
                K.steps = 1;
                K.scheme = schemeName;
               
            elseif ( strcmpi(schemeName, 'minchev') == 1 )
                K.order = 4;
                K.internalStages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'strehmel-weiner') == 1 )
                K.order = 4;
                K.internalStages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % EPI RUNGE-KUTTA:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LAWSON:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'lawson4') == 1 )
                K.order = 4;
                K.internalStages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % GENERALIZED LAWSON:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'genlawson43') == 1 )
                K.order = 4;
                K.internalStages = 4;
                K.steps = 3;
                K.scheme = schemeName;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MODIFIED GENERALIZED LAWSON:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'modgenlawson43') == 1 )
                K.order = 4;
                K.internalStages = 4;
                K.steps = 3;
                K.scheme = schemeName;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PREDICTOR-CORRECTOR:    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'pecec433') == 1 )
                K.order = 4;
                K.internalStages = 3;
                K.steps = 3;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'pecec635') == 1 )
                K.order = 6;
                K.internalStages = 3;
                K.steps = 5;
                K.scheme = schemeName;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MISC:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'eglm433') == 1 )
                K.order = 4;
                K.internalStages = 3;
                K.steps = 3;
                K.scheme = schemeName;
                
            else
                error('SPINSCHEME:constructor', 'Unrecognized scheme.')
            end
            
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE AND STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = true )
        
        %  Evaluate a phi function:
        phi = phiEval(l, LR, N, dim, nVars)
        
        % Evaluate a phit function:
        phit = phitEval(l, C, LR, N, dim, nVars)
        
        % Get a function handle to a phi function:
        phi = phiFun(l)
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE AND NON-STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = false )
       
        % Compute coefficients of a SPINSCHEME:
        coeffs = computeCoeffs(K, dt, L, LR, S)
        
        % Do one step of a SPINSCHEME:
        uSol = oneStep(K, coeffs, Nc, S, uSol)
        
        % Get enough initial data when using a multistep SPINSCHEME:
        [uSol, dt] = startMultistep(K, adapTime, dt, L, LR, Nc, pref, S, uSol)
    end
    
end