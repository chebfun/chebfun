classdef coeffsDiscretization < opDiscretization
%COEFFSDISCRETIZATION   Abstract class for coefficient-based discretization of 
%operators.
%   COEFFSDISCRETIZATION is a partial implementation of OPDISCRETIZATION using
%   Chebyshev coefficients or Fourier coefficients. COEFFSDISCRETIZATION cannot 
%   be used directly as a discretization for LINOPs. ULTRAS and TRIGSPEC are 
%   full implementations.
%
% See also OPDISCRETIZATION, TRIGSPEC, ULTRAS.

%  Copyright 2016 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = coeffsDiscretization(source, dimension, dom)
            
            if ( (nargin == 0) || isempty(source) )
                % Construct an empty COEFFSDISCRETIZATION.
                return
            end
            
            % Attach SOURCE and the DOMAIN information to the object:
            disc.source = source;
            disc.domain = source.domain;
            
            % Obtain the coeffs and output space required for this source:
            disc.coeffs = disc.getCoeffs(source);
            
            % Determine the dimension adjustments and outputSpace:
            disc.dimAdjust = disc.getDimAdjust(source);
            disc.projOrder = disc.getProjOrder(source);
            disc.outputSpace = disc.getOutputSpace(source);
            
            % Assign DIMENSIONS and DOMAIN if they were passed.
            if ( nargin > 1 )
                disc.dimension = dimension;
            end
            if ( nargin > 2 )
                disc.domain = domain.merge(dom, disc.domain);
            end
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true)
        
        function dimVals = dimensionValues(pref)
            %DIMENSIONVALUES    Return a vector of desired discretization sizes.
            %  DIMVALS = DIMENSIONVALUES(PREF) returns a vector containing
            %  elements that prescribe what values should be used as dimension
            %  values for discretizating linear operators. DIMVALS is affected
            %  by the minimum and maximum discretizations specified in the
            %  CHEBOPPREF object PREF.
            
            % We simply go up in powers of 2.
            
            minPow = log2(pref.minDimension);
            maxPow = log2(pref.maxDimension);
            
            if ( minPow > maxPow )
                error('CHEBFUN:COEFFSDISCRETIZATION:dimensionValues', ...
                    ['Minimum discretiation specified is greater than ' ...
                    'maximum discretization specified']);
            end
            
            % Go up in powers of 2:
            powVec = minPow:maxPow;
            dimVals = round(2.^powVec);
            
        end
        
        % Get coefficient representation of the source.
        c = getCoeffs(source)
        
        % Obtain the range of the spectral operator.
        outputSpace = getOutputSpace(source)
        
        % Sparse Toeplitz matrix.
        T = sptoeplitz(col, row)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS IMPLEMENTED IN THIS FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Access = public, Static = false )
        function varargout = convertOperator(varargin)
            error('COEFFSDISCRETIZATION:instantiate:fail', ...
                ['Cannot represent this operator. ', ...
                 'Suggest you use VALSDISCRETIZATION.'])
        end
    end

end
