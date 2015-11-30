classdef ultraS < coeffsDiscretization
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

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
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
        
        function disc = ultraS(varargin)
            disc = disc@coeffsDiscretization(varargin{:});
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
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true)
        
        % Compute sparse representation for conversion operators.
        T = spconvert(n, lam)
        
        % Construct a sparse Hankel operator.
        H = sphankel(r)
        
    end
    
end
