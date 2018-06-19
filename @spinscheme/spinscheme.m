classdef spinscheme
%SPINSCHEME   Abstract class for representing time-stepping schemes.
%   SPINSCHEME is an abstract class for representing time-stepping schemes.
%   Full implementations are EXPINTEG (exponential integrators) and IMEX 
%   (implicit-explcit schemes).
%
% See also EXPINTEG and IMEX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        
        stages          % number of internal stages (1x1 INT)
        order           % order of the method (1x1 INT)
        scheme          % Time-stepping scheme (STRING)
        steps           % number of previous time-steps used, 1 if one-step
                        % method, > 1 if multistep method (1x1 INT)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT AND NON-STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = false )
        
        % Compute coefficients of a SPINSCHEME:
        schemeCoeffs = computeCoeffs(K, dt, L, M, S)
        
        % Do one step of a SPINSCHEME:
        [uSol, NuSol] = oneStep(K, dt, schemeCoeffs, Nc, Nv, nVars, S, uSol, NuSol)
        
        % Get enough initial data when using a multistep SPINSCHEME:
        [uSol, NuSol] = startMultistep(K, dt, L, Nc, Nv, pref, S, uSol, NuSol)
        
    end
    
end
