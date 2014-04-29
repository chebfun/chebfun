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
    
%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    properties
        % Coefficients of the operator
        coeffs
        % The range of the ultrapspherical spectral operator.
        outputSpace = [];
    end
    
    methods
        function disc = ultraS(source, dimension, domain)
            %ULTRAS(SOURCE, DIMENSION, DOMAIN)   ULTRAS constructor.
            
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
            disc.domain = chebfun.mergeDomains(source.domain, disc.domain); 
            
            % Obtain the coeffs and output space required for this source:
            disc.coeffs = ultraS.getCoeffs(source);
            disc.outputSpace = ultraS.getOutputSpace(source);
            
        end
        
        
%         function [isDone, epsLevel] = testConvergence(disc,v)
%             % Test convergence on a single interval:
%             v = full(v);
%             f = chebtech2({[], flipud(v)});
%             [isDone, epsLevel] = strictCheck(f);
%         end

    end
    
    methods ( Access = private )
        
        % Conversion (transformation) operator for Ultraspherical method.
        S = convert(A, K1, K2)

    end
    
    methods ( Access = private, Static = true)
        
        % Conversion matrix used in the ultraspherical spectral method.
        S = convertmat(n, K1, K2)
        
        % Differentiation matrices for ultraspherical spectral method.
        D = diffmat(n, m)
        
        % Get coefficients.
        c = getCoeffs( source )
        
        % Obtain the range of the ultrapspherical spectral operator.
        outputSpace = getOutputSpace(source)
        
        % Compute sparse representation for conversion operators. 
        T = spconvert(n, lam)
        
        % Construct a sparse Hankel operator.
        H = sphankel(r)

        % Sparse Toeplitz matrix.
        T = sptoeplitz(col, row)

    end
end
