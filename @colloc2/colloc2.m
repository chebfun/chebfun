classdef colloc2 < chebDiscretization
    
    properties (Access=private)
        mldivideData = [];
    end
     
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    methods
        function disc = colloc2(source,dimension,domain)
            if isempty(source)
                return
            end
            disc.source = source; 
            disc.domain = source.domain;

            if nargin > 1
                disc.dimension = dimension;
                if nargin > 2
                    disc.domain = domain;
                end
            end
        end
         
    end
        
end