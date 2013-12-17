classdef chebfun2
    
    properties ( Access = public )
        cols
        rows
        pivotValues
        domain = [-1 1 -1 1];
    end
    
    methods
        
        function f = chebfun2(varargin)
            % The main CHEBFUN2 constructor!
            
            % Return an empty CHEBFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            f = constructor(f, varargin{:});
            
        end
        
    end
    
    % Static methods implemented by CHEBFUN class.
    methods ( Static = true )
        
        X = coeffs2vals( U ); 
        
        X = vals2coeffs( U ); 
        
        [xx, yy] = chebpts2(nx, ny, domain);
        
        F = outerProduct(f, g);   % outer-product of two chebfuns.
    end

    % Private methods implemented by CHEBFUN2 class.
    methods ( Access = private )
        
    end
    
    % Static private methods implemented by CHEBFUN2 class.
    methods ( Static = true, Access = private )
        
    end
    
    % Methods implemented by CHEBFUN2 class.
    methods
        
    end
    
end

