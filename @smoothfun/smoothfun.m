classdef smoothfun < onefun % (Abstract) 
%SMOOTHFUN   Approximate smooth functions on [-1,1]. 
%
%   Abstract (interface) class for approximating smooth functions on the
%   interval [-1,1].
%
% Constructor inputs:
%   SMOOTHFUN.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF) constructs a SMOOTHFUN
%   object on the interval [-1,1] from the function handle OP. Currently the
%   only subclass of SMOOTHFUN is CHEBTECH, so SMOOTHFUN will call
%   CHEBTECH.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF2), where PREF2 is PREF merged
%   with the default CHEBTECH preferences.
%
% See also SMOOTHFUN.PREF, CHEBTECH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHFUN Class Description:
%
% The SMOOTHFUN class is an abstract class for representations of smooth
% functions on the interval [-1,1].
%
% Currently the only types of SMOOTHFUNs are CHEBTECH objects, which represent
% the smooth functions by Chebyshev interpolants.
%
% Class diagram: [<<onefun>>] <-- [<<SMOOTHFUN>>] <-- [<<chebtech>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Constructor for the SMOOTHFUN class.
    methods (Static)
        function obj = constructor(op, vscale, hscale, pref)
            
            % We can't return an empty SMOOTHFUN, so pass an empty OP down.
            if ( nargin == 0 )
                op = [];
            end
            

            % Define vscale if none given:
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            % Define hscale if none given:
            if ( nargin < 3 || isempty(hscale) )
                hscale = 1;
            end
            % Determine preferences if not given, merge if some are given:
            if ( nargin < 4 || isempty(pref) )
                pref = smoothfun.pref;
            else
                pref = smoothfun.pref(pref);
            end
            
            % Merge preferences:
            pref = chebtech.pref(pref, pref.smoothfun);
            % Call the CHEBTECH constructor
            obj = chebtech.constructor(op, vscale, hscale, pref);
            
        end
        
    end
    
    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY SMOOTHFUN CLASS.
    methods ( Abstract = true )

    end

    %% ABSTRACT STATIC METHODS REQUIRED BY SMOOTHFUN CLASS.
    methods ( Abstract = true, Static = true )
        
    end
    
    %% Methods implimented by SMOOTHFUN class.
    methods 
        
    end
    
    %% Static methods implimented by SMOOTHFUN class.
    methods ( Static = true ) 
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
    end
    
end


