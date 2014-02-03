classdef ultraS < chebDiscretization
%ULTRAS   ULTRASPHERICAL class for discretizating differential operators 
%         on bounded domain. 
%
%    This class uses the ultraspherical spectral to discretize and ultimately
%    solve 1D boundary value problem defined on a bounded interval.  
% 
% ULTRAS(SOURCE, DIMENSION, DOMAIN) constructs a ultraspherical spectral 
% discretization of SOURCE of size DIMENSION-by-DIMENSION for BVPs on the
% interval DOMAIN. 
% 
% ULTRAS(SOURCE, DIMESION) takes the DOMAIN from SOURCE. 
%
% ULTRAS(SOURCE) takes the dimension from SOURCE. 

% For more details about the ultrapsherical spectral methods, see: 
% S. Olver and A. Townsend, A fast and well-conditioned spectral method, SIAM
% Review, 55 (2013), pp. 462-489.
    
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
            disc.domain = chebfun.mergeDomains(source.domain,disc.domain); 
            
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
