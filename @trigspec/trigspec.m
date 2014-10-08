classdef trigspec < chebDiscretization
%TRIGSPEC    Fourier spectral method in coefficient space.
%   TRIGSPEC is an implementation of CHEBDISCRETIZATION that implements a
%   Fourier spectral method in coefficient space.
%
% See also TRIGCOLLOC.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        coeffs        % Coefficients of the operator.
        outputSpace   % The range of the operator.
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = trigspec(source, dimension, dom)
            
            if ( (nargin == 0) || isempty(source) )
                % Construct an empty TRIGSPEC object.
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
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        function tech = returnTech()
            %RETURNTECH    Return the appropriate tech to use for TRIGSPEC.
            tech = @fourtech;
        end
        
        % Differentiation matrices for TRIGSPEC.
        D = diffmat(n, m)
        
        % Multiplication matrices for TRIGSPEC.
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
    methods ( Access = private, Static = true )
        
        % Get coefficient representation of the source.
        c = getCoeffs(source)
        
        % Obtain the range of the operator.
        outputSpace = getOutputSpace(source)
     
        % Construct a sparse Hankel operator.
        H = sphankel(r)
        
        % Sparse Toeplitz matrix.
        T = sptoeplitz(col, row)
        
    end
    
end
