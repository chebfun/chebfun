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
%
% For more details about the ultrapsherical spectral methods, see: S. Olver and
% A. Townsend, A fast and well-conditioned spectral method, SIAM Review, 55
% (2013), pp. 462-489.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        coeffs        % Coefficients of the operator
        outputSpace   % The range of the ultraspherical spectral operator
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = ultraS(source, dimension, dom)
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
                disc.domain = domain.merge(dom, disc.domain);
            end
            
        end
        
        % Dimension reduction for operator matrix.
        [PA, P, PS] = reduce(disc, A, S)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE METHODS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = false )
        
        % Conversion (transformation) operator for Ultraspherical method.
        S = convert(A, K1, K2)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true)
        
        function tech = returnTech()
            %RETURNTECH    Return the appropriate tech to use for ULTRAS.
            tech = @chebtech2;
        end
        
        % Conversion matrix used in the ultraspherical spectral method.
        S = convertmat(n, K1, K2)
        
        % Differentiation matrices for ultraspherical spectral method.
        D = diffmat(n, m)
        
        % Multiplication matrices for ultraspherical spectral method.
        D = multmat(n, f, lambda)
        
        function dimVals = dimensionValues(pref)
            %DIMENSIONVALUES   Return a vector of desired discretization sizes.
            %  DIMVALS = DIMENSIONVALUES(PREF) returns a vector containing
            %  elements that prescribe what values should be used as dimension
            %  values for discretizating linear operators. DIMVALS is affected
            %  by the minimum and maximum discretizations specified in the
            %  CHEBOPPREF object PREF.
            
            % We simply go up in powers of 2.
            
            minPow = log2(pref.minDimension);
            maxPow = log2(pref.maxDimension);
            
            if ( minPow > maxPow )
                error('CHEBFUN:ULTRAS:ultraS:dimensionValues', ...
                    ['Minimum discretiation specified is greater than ' ...
                    'maximum discretization specified']);
            end
            
            % Go up in powers of 2:
            powVec = minPow:maxPow;
            dimVals = round(2.^powVec);
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true)
        
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
