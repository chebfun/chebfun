classdef abstractfun % (Abstract) 
%ABSTRACTFUN   Approximate global functions.
% Class diagram: [<<Chebfun>>] <>-- [<<ABSTRACTFUN>>] <-- [<<fun>>]
%                                                     <-- [deltafun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    methods (Static)
        function obj = constructor(op, vscale, hscale, pref)
            
            % We can't return an empty ABSTRACTFUN, so pass an empty OP down.
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
                pref = chebpref();
            else
                pref = chebpref(pref);
            end

            % Call the relevent constructor:
            if ( isa(op, 'deltafun') )
                % OP is already a DELTAFUN!
                obj = op;
                
            elseif ( pref.deltaEnabled )
                % Delta function mode, call DELTAFUN constructor:
                obj = deltafun(op, [], [], vscale, hscale, pref);

                % Return just a SMOOTHFUN if no singularities found:
                if ( ~hasdeltas(obj) )
                    obj = obj.funPart; 
                end 
                
            else
                % STANDARD mode; call FUN constructor:
                obj = fun.constructor(op, vscale, hscale, pref);

            end
        
        end
        
    end
    
    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY ABSTRACTFUN CLASS.
    methods ( Abstract = true )
        
            end

    %% ABSTRACT STATIC METHODS REQUIRED BY ABSTRACTFUN CLASS.
    methods ( Abstract = true, Static = true )

    end
    
    %% Methods implemented by ABSTRACTFUN class.
    methods 
        
    end
    
    %% Static methods implemented by ABSTRACTFUN class.
    methods ( Static = true ) 
        
    end
    
end
