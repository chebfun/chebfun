classdef colloc2 < colloc
    
    methods
        function disc = colloc2(varargin)
 %COLLOC2    Collocation discretization on 2nd kind points.
            
            disc = disc@colloc(varargin{:});
       end
    end
    
end