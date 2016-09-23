classdef opDiscretization 
%OPDISCRETIZATION    Convert a chebmatrix or linop to discrete form.
%   This class is not called directly by the end user. 
%
% See also COLLOC, ULTRAS.

%  Copyright 2016 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%   Objects of this class store a source (either a CHEBMATRIX or a LINOP),
%   domain, and dimension (vector of discretization lengths). In other words,
%   they have a CHEBMATRIX or LINOP object, which represent an abstract linear
%   operator (with or without constraints). Most importantly, objects of types
%   of concrete implementation of OPDISCRETIZATION are able to instantianate
%   themselves to a matrix corresponding to discretization of the operator on a
%   Chebyshev grid. Objects of this type also have to implement a mldivide
%   (backslash) method, which yields solution to problems of ODEs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        source = []       % linop or chebmatrix to be discretized
        domain = []       % may generalize that of the source
        dimension = []    % vector of lengths, one per subinterval
        dimAdjust = []    % size of the input space relative to disc.dimension
        projOrder = []    % projection order for rectangularizing
    end
        
    properties ( Dependent = true )
        numIntervals      % number of intervals in the domain
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS: (ABSTRACT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = true, Static = false )   
        
        % Converts a CHEBFUN into a vector of values (or coefficients,
        % depending on the implementation). 
        values = toValues(disc, f)
        
        % Converts a vector of values (or coefficients) to a CHEBFUN.
        f = toFunctionIn(disc, values)
        
        % Converts a vector of values (or coefficients) to a CHEBFUN.
        f = toFunctionOut(disc, values)
        
        % Returns a linear system RHS using the designated discretization
        % parameters.
        b = rhs(disc, f, varargin)
        
        % Reduces (projects) block rows to make space for the constraints.
        [PA, P, PS] = reduce(disc, blocks)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS: (ABSTRACT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = true, Static = true )
        
        % Return a vector of desired discretization sizes.
        dimVals = dimensionValues(pref)
        
        % Return the appropriate tech to be used with the discretization.
        tech = returnTech()
        
        % Get dimension adjustment for EXPM.
        expmDimAdjust = getExpmDimAdjust(L)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS: (CONCRETE)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Get dimension adjustment:
        dimAdjust = getDimAdjust(L)
        
        % Get projection order:
        projOrder = getProjOrder(L)
        
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS IMPLEMENTED IN THIS FILE:    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        function n = get.numIntervals(disc)
            %NUMINTERVALS   Number of subintervals a OPDISCRETIZATION acts on.
            n = length(disc.domain) - 1;
        end
        
        function t = isempty(disc)
            %ISEMPTY   Check if source property of a OPDISCRETIZATION is empty.
            t = isempty(disc.source);
        end
        
        function t = isFactored(disc) %#ok<MANU>
            %ISFACTORED   Check if a factorization of the source already exists.
            %   This method gives a discretization a chance to overload and store
            %   matrix factors for the purpose of short-circuiting the linsolve
            %   process. By default it never happens.
            t = false;
        end
        
        function [x, disc] = mldivide(disc, A, b)
            %OPDISCRETIZATION.MLDIVIDE
            %   By default, the solution of a discrete Ax = b uses standard
            %   backslash. But concrete implementations may overload it.
            x = A\b;
        end
        
    end        
    
end
