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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    properties
        coeffs        % Coefficients of the operator
        outputSpace   % The range of the ultraspherical spectral operator
    end
    
    methods
        function disc = ultraS(source, dimension, domain)
            %ULTRAS(SOURCE, DIMENSION, DOMAIN)   ULTRAS constructor.
            
            if ( (nargin == 0) || isempty(source) )
                % Construct an empty ULTRAS.
                return
            end
            
            % Attach SOURCE and the DOMAIN information to the object:
            disc.source = source; 
            disc.domain = source.domain;
            % Obtain the coeffs and output space required for this source:
            disc.coeffs = ultraS.getCoeffs(source);
            % Determine the dimension adjustments and outputSpace:
            disc.dimAdjust = ultraS.getDimAdjust(source);
            disc.projOrder = ultraS.getProjOrder(source);
            disc.outputSpace = ultraS.getOutputSpace(source);
            
            % Assign DIMENSIONS and DOMAIN if they were passed.
            if ( nargin > 1 )
                disc.dimension = dimension;
            end
            if ( nargin > 2 )
                disc.domain = chebfun.mergeDomains(domain, disc.domain); 
            end
            
        end

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
        
        % Get coefficient representation of the source.
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
