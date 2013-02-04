classdef (Abstract) smoothfun < onefun
    
    % Copyright 2011 by The University of Oxford and The Chebfun Developers.
    % See http://www.chebfun.org for Chebfun information. 
    
    methods (Static)
        
        function obj = constructor(op, pref)
            % Constructor for the SMOOTHFUN class.
            
            if nargin == 0
                return
            end
            
            % Call the relevent constructor
            if nargin < 2
                pref = smoothfun.pref;
            else
                pref = smoothfun.pref(pref);
            end

            if ( pref.smoothfun.chebkind == 1 )
                pref = fun1.pref(pref, pref.smoothfun);
                obj = fun1(op, [], [], pref);
            else
                pref = fun2.pref(pref, pref.smoothfun);
            
                obj = fun2(op, [], [], pref);
            end
            
        end
        
    end
    
    methods % Methods implimented by SMOOTHFUN class.
        
    end
    
    methods ( Static = true ) % Static methods implimented by SMOOTHFUN class.
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
    end
    
end


