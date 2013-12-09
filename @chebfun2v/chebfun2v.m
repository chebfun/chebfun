classdef chebfun2v
    
    properties ( Access = public )
        components   % Array of chebfun2 objects. 
        nComponents  % Number of components
    end
    
    methods
        
        function f = chebfun2v(varargin)
            % The main CHEBFUN2V constructor!
            
            % Return an empty CHEBFUN2V:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            f = constructor(f, varargin{:});
            
        end
    end 
    
end