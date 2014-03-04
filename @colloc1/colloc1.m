classdef colloc1 < colloc
    
    methods
        function disc = colloc1(varargin)
 %COLLOC1    Collocation discretization on 1st kind points.
            
            disc = disc@colloc(varargin{:});
       end
    end
    
end