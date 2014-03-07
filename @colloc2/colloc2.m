classdef colloc2 < colloc
%   operators and systems. 
%   The default discretization type to use is set by CHEBOPPREF. You can also
%   use CHEBOPPREF to create a preferences object and change its
%   'discretization' property. 
%   flops, where in most contexts N is determined automatically to resolve the
%   solution. The allowed values of N are governed by the 'dimensionValues'
%   property in CHEBOPPREF. You can also set the maximum N (including systems
%   and piecewise definitions) through the 'maxTotalLength' property.
%   See also CHEBOPPREF, CHEBDISCRETIZATION.
    
    methods
        function disc = colloc2(varargin)
 %COLLOC2    Collocation discretization on 2nd kind points.
            
            disc = disc@colloc(varargin{:});
       end
    end
    
end