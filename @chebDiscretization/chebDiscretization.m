classdef (Abstract) chebDiscretization 
%CHEBDISCRETIZATION Convert a chebmatrix or linop to discrete form.
%   This class is not called directly by the end user. 
%
%   See also COLLOC, ULTRAS.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes:
%
% Objects of this class store a source (either a CHEBMATRIX or a LINOP), domain,
% and dimension (vector of discretization lengths). In other words, they have a
% CHEBMATRIX or LINOP object, which represent an abstract linear operator (with
% or without constraints). Most importantly, objects of types of concrete
% implementation of CHEBDISCRETIZATION are able to instantianate themselves to a
% matrix corresponding to discretization of the operator on a Chebyshev grid.
% Objects of this type also have to implement a mldivide (backslash) method,
% which yields solution to problems of ODEs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    properties
        source = []       % linop or chebmatrix to be discretized
        domain = []       % may generalize that of the source
        dimension = []    % vector of lengths, one per subinterval
        dimAdjust = []    % size of the input space relative to disc.dimension
    end
        
    properties ( Dependent )
        numIntervals      % number of intervals in the domain
    end

    methods
        
        function n = get.numIntervals(disc)
        %NUMINTERVALS    Number of subintervals a CHEBDISCRETIZATION acts on.
            n = length(disc.domain) - 1;
        end        


        function t = isempty(disc)
        % Returns true if source property of a CHEBDISCRETIZATION is empty.
            t = isempty(disc.source);
        end

        function t = isFactored(disc)
        %CHEBDISCRETIZATION.ISFACTORED
        %
        % This method gives a discretization a chance to overload and store
        % matrix factors for the purpose of short-circuiting the linsolve
        % process. By default it never happens.
            t = false;
        end
        
        function [x, disc] = mldivide(disc, A, b)
        %CHEBDISCRETIZATION.MLDIVIDE 
        %
        % By default, the solution of a discrete Ax=b uses standard backslash.
        % But concrete implementations may overload it.
            x = A\b;
        end           
        
    end
    
    methods ( Static )
        space = getDimAdjust(L)
    end
            
    methods ( Abstract )        
        % Converts a chebfun into a vector of values (or coefficients,
        % depending on the implementation). 
        values = toValues(disc, f)
        
        % Converts a vector of values (or coefficients) to a chebfun.
        f = toFunction(disc, values)
        
        % Returns a linear system RHS using the designated discretization
        % parameters.
        b = rhs(disc, f, varargin)
        
        % Reduces (projects) block rows to make space for the constraints.
        [PA, P, PS] = reduce(disc, blocks)
        
        % Extract the j-k block from a discretization.
        discjk = extractBlock(disc, j, k)
        
    end
    
end
