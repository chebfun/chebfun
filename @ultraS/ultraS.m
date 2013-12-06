classdef ultraS < chebDiscretization
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    properties
        coeffs
        outputSpace = [];
    end
    
    methods
        function disc = ultraS(source,dimension,domain)
            %ULTRAS constructor.
            
            if ( nargin == 0 || isempty(source) )
                % Construct an empty ULTRAS.
                return
            end
            
            if ( nargin > 1 )
                disc.dimension = dimension;
                if ( nargin > 2 )
                    disc.domain = domain;
                end
            end
            
            disc.source = source; 
            disc.domain = source.domain;
            
            % Obtain the coeffs and output psace required for this source:
            disc.coeffs = getCoeffs(source);
            disc.outputSpace = getOutputSpace(source);
            
        end
        
        
%         function [isDone, epsLevel] = testConvergence(disc,v)
%             % Test convergence on a single interval:
%             v = full(v);
%             f = chebtech2({[], flipud(v)});
%             [isDone, epsLevel] = strictCheck(f);
%         end

    end
    
end
