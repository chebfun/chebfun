classdef expinteg < spinscheme
%EXPINTEG   Class for representing exponential integrators.
%   EXPINTEG is a class for representing exponential integrators, which are used 
%   for the time integration in SPIN, SPIN2 and SPIN3.
%
% Construction: 
%
%   K = EXPINTEG(SCHEME) constructs an EXPINTEG object corresponding to the 
%   exponential integrator SCHEME. SCHEME is a (case-insensitive) STRING.
%
% Available exponential integrators:
%
%   ETD ADAMS-BASHFORTH: 'abnorsett4', 'abnorsett5', 'abnorsett6'
%
%   ETD RUNGE-KUTTA: 'etdrk4', 'exprk5s8', 'friedli', 'hochbruck-ostermann',
%                    'krogstad', 'minchev', 'strehmel-weiner'
%
%   LAWSON: 'lawson4', 'ablawson4'
%
%   GENERALIZED LAWSON: 'genlawson41', genlawson42', genlawson43', genlawson44',
%                       'genlawson45'
%
%   MODIFIED GENERALIZED LAWSON: 'modgenlawson41', modgenlawson42', 
%                                'modgenlawson43', modgenlawson44',
%                                'modgenlawson45'
%
%   PREDICTOR-CORRECTOR: 'pec423', 'pecec433', 'pec524', 'pecec534', 'pec625',
%                        'pecec635', 'pec726', 'pecec736'
%
% See also SPIN, SPIN2, SPIN3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function K = expinteg(schemeName)
            
            if ( nargin == 0 )
                return
            end
                           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ETD ADAMS-BASHFORTH:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ( strcmpi(schemeName, 'abnorsett4') == 1 )
                K.order = 4;
                K.stages = 1;
                K.steps = 4;
                K.scheme = schemeName;

            elseif ( strcmpi(schemeName, 'abnorsett5') == 1 )
                K.order = 5;
                K.stages = 1;
                K.steps = 5;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'abnorsett6') == 1 )
                K.order = 6;
                K.stages = 1;
                K.steps = 6;
                K.scheme = schemeName;
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ETD RUNGE-KUTTA:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'etdrk2') == 1 )
                K.order = 2;
                K.stages = 2;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'etdrk4') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'exprk5s8') == 1 )
                K.order = 5;
                K.stages = 8;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'friedli') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'hochbruck-ostermann') == 1 )
                K.order = 4;
                K.stages = 5;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'krogstad') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'minchev') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'strehmel-weiner') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LAWSON:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'ablawson4') == 1 )
                K.order = 4;
                K.stages = 1;
                K.steps = 4;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'lawson4') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % GENERALIZED LAWSON:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'genlawson41') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'genlawson42') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 2;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'genlawson43') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 3;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'genlawson44') == 1 )
                K.order = 5;
                K.stages = 4;
                K.steps = 4;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'genlawson45') == 1 )
                K.order = 6;
                K.stages = 4;
                K.steps = 5;
                K.scheme = schemeName;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MODIFIED GENERALIZED LAWSON:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'modgenlawson41') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 1;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'modgenlawson42') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 2;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'modgenlawson43') == 1 )
                K.order = 4;
                K.stages = 4;
                K.steps = 3;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'modgenlawson44') == 1 )
                K.order = 5;
                K.stages = 4;
                K.steps = 4;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'modgenlawson45') == 1 )
                K.order = 6;
                K.stages = 4;
                K.steps = 5;
                K.scheme = schemeName;
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PREDICTOR-CORRECTOR:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ( strcmpi(schemeName, 'pec423') == 1 )
                K.order = 4;
                K.stages = 2;
                K.steps = 3;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'pecec433') == 1 )
                K.order = 4;
                K.stages = 3;
                K.steps = 3;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'pec524') == 1 )
                K.order = 5;
                K.stages = 2;
                K.steps = 4;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'pecec534') == 1 )
                K.order = 5;
                K.stages = 3;
                K.steps = 4;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'pec625') == 1 )
                K.order = 6;
                K.stages = 2;
                K.steps = 5;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'pecec635') == 1 )
                K.order = 6;
                K.stages = 3;
                K.steps = 5;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'pec726') == 1 )
                K.order = 7;
                K.stages = 2;
                K.steps = 6;
                K.scheme = schemeName;
                
            elseif ( strcmpi(schemeName, 'pecec736') == 1 )
                K.order = 7;
                K.stages = 3;
                K.steps = 6;
                K.scheme = schemeName;
  
            else
                error('EXPINTEG:constructor', 'Unrecognized scheme.')
            end
        
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE AND STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = true )
        
        % Get a function handle to a gamma-function:
        g = gammaFun(j, k)
        
        %  Evaluate a phi-function:
        g = gammaEval(j, k, LR, N, dim, nVars)
        
        %  Evaluate a phi-function:
        phi = phiEval(l, LR, N, dim, nVars)
        
        % Get a function handle to a phi-function:
        phi = phiFun(l)
        
        % Evaluate a psi-function:
        psi = psiEval(l, C, LR, N, dim, nVars)
        
    end
    
end
