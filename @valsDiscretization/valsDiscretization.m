classdef valsDiscretization < opDiscretization
%VALSDISCRETIZATION   Abstract class for collocation discretization of 
%operators.
%   VALSDISCRETIZATION is a partial implementation of OPDISCRETIZATION using
%   first-kind Chebyshev points, second-kind Chebyshev points, or equally
%   spaced points. VALSDISCRETIZATION cannot be used directly as a 
%   discretization for LINOPs. CHEBCOLLOC1, CHEBCOLLOC2, and TRIGCOLLOC are full 
%   implementations.
%
% See also CHEBCOLLOC, TRIGCOLLOC, CHEBCOLLOC1, CHEBCOLLOC2, OPDISCRETIZATION.

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = private )
        % Stores LU factors of a matrix, for repeated solves at fixed size:
        mldivideData = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = valsDiscretization(source, dimension, domain)
            %VALSDISCRETIZATION    Collocation discretization constructor.
            %   (Called by subclasses for parts in common.)
            
            if ( (nargin == 0) || isempty(source) )
                % Construct an empty COLLOC.
                return
            end
            
            % Attach SOURCE and the DOMAIN information to the object:
            disc.source = source;
            disc.domain = source.domain;
            % Determine the dimension adjustment:
            disc.dimAdjust = valsDiscretization.getDimAdjust(source);
            disc.projOrder = valsDiscretization.getProjOrder(source);
            
            % Assign DIMENSIONS and DOMAIN if they were passed:
            if ( nargin > 1 )
                disc.dimension = dimension;
            end
            if ( nargin > 2 )
                disc.domain = domain;
            end
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Abstract = true )

        % Indefinite integration:
        C = cumsum(disc)
        
        % Differentiation:
        D = diff(disc, m)
        
        % Points where function values are represented:
        [x, w, v, t] = functionPoints(disc)
        
        % Points where equations are enforced:
        [x, w, v, t] = equationPoints(disc)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS (IMPLEMENTED BY THIS ABSTRACT CLASS.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Discretization points: used by CHEBCOLLOC1, CHEBCOLLOC2, and 
        % TRIGCOLLOC.
        [x, w, v, t] = points(varargin);
        
        function dimVals = dimensionValues(pref)
            %DIMENSIONVALUES   Return a vector of desired discretization sizes.
            %  DIMVALS = DIMENSIONVALUES(PREF) returns a vector containing
            %  elements that prescribe what values should be used as dimension
            %  values for discretizating linear operators. DIMVALS is affected
            %  by the minimum and maximum discretizations specified in the
            %  CHEBOPPREF object PREF.
            
            % We want to go up in powers of 2 up to a point, after which, we go
            % up in half powers of two.
            minPow = log2(pref.minDimension);
            maxPow = log2(pref.maxDimension);
            
            if ( minPow > maxPow )
                error('CHEBFUN:VALSDISCRETIZATION:valsDiscretization:dimensionValues', ...
                    ['Minimum discretiation specified is greater than ' ...
                     'maximum discretization specified.']);
            end
            
            if ( maxPow <= 9 )
                % We're happy to go up in steps of 2 up until 512.
                powVec = minPow:maxPow;
            elseif ( minPow >= 9 )
                powVec = minPow:.5:maxPow;
            else
                powVec = [minPow:9, 9.5:.5:maxPow];
            end
            
            dimVals = round(2.^powVec);
            
        end
        
    end
    
end
