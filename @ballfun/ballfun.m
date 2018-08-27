classdef ballfun
    %BALLFUN  A library for the DFS method in the solid sphere.
    
    properties
        
        coeffs   % Chebyshev-Fourier-Fourier coefficients array of a BALLFUN function
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = ballfun(varargin)
            % The main BALLFUN constructor!
            
            % Return an empty BALLFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            % Call the constructor, all the work is done here:
            f = constructor(f, varargin{:});
       
        end
    end
    
    methods ( Access = public, Static = true )
        
        % Convert to Chebyshev--Fourier--Fourier values
        VALS = coeffs2vals(CFS);
        
        % Convert to Chebyshev--Fourier--Fourier values
        CFS = vals2coeffs(VALS);
        
        % Convert a chebfun function to a ballfun function
        f = chebfun2ballfun(g,S);
        
        % Convert a spherefun function to a ballfun function
        f = spherefun2ballfun(g,S);
        
        % Compute the solid harmonics
        f = solHarm(l,m);
        
    end
    
    methods ( Access = private, Static = true )
        
    end
    
    methods ( Access = private, Static = false )
        
    end
    
end
