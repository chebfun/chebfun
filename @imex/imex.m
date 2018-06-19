classdef imex < spinscheme
%IMEX   Class for representing implicit-explicit schemes.
%   IMEX is a class for representing implicit-explicit schemes, which are used 
%   for the time integration in SPINSPHERE.
%
% Construction: 
%
%   K = IMEX(SCHEME) constructs an IMEX object corresponding to the 
%   IMEX scheme SCHEME. SCHEME is a (case-insensitive) STRING.
%
% Available IMEX schemes:
%
%   IMEX BDF: 'imexbdf4' (for diffusive PDEs only)
%
%   IMEX RUNGE-KUTTA: 'lirk4' (for both diffusive and dispersive PDEs)
%
% See also SPINSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function K = imex(schemeName)
            
            if ( nargin == 0 )
                return
            end
   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % IMEX RUNGE-KUTTA:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ( strcmpi(schemeName, 'imexbdf4') == 1 )
                K.order = 4;
                K.stages = 1;
                K.steps = 4;
                K.scheme = schemeName;
            
            elseif ( strcmpi(schemeName, 'lirk4') == 1 )
                K.order = 4;
                K.stages = 5;
                K.steps = 1;
                K.scheme = schemeName;
                
            else
                error('IMEX:constructor', 'Unrecognized scheme.')
            end
        
        end
        
    end
    
end
