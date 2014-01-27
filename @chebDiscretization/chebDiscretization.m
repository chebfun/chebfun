classdef (Abstract) chebDiscretization 
%CHEBDISCRETIZATION Convert a chebmatrix or linop to discrete form.
%   This class is not called directly by the end user. 
%
%   See also COLLOC2, ULTRAS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.
    
% Objects of this class store a source (either a chebmatrix or a linop),
% domain, and dimension (vector of discretization lengths). 
    
    properties
        source = []       % linop or chebmatrix to be discretized
        domain = []       % may generalize that of the source
        dimension = []    % vector of lengths, one per subinterval
    end
        
    properties (Dependent)
        numIntervals      % number of intervals in the domain
    end

    methods
        function n = get.numIntervals(disc)
            n = length(disc.domain) - 1;
        end        

        function t = isempty(disc)
            t = isempty(disc.source);
        end

        % This method gives a discretization a chance to overload and store
        % matrix factors for the purpose of short-circuiting the linsolve
        % process. By default it never happens.
        function t = isFactored(disc)
            t = false;
        end
        
        % By default, the solution of a discrete Ax=b uses standard
        % backslash. But implementations may overload it. 
        function [x,disc] = mldivide(disc,A,b)
            x = A\b;
        end           
        
    end
        
    methods (Abstract)        
        % Converts a chebfun into a vector of values (or coefficients,
        % depending on the implementation). 
        values = toValues(disc,f)
        
        % Converts a vector of values (or coefficients) to a chebfun.
        f = toFunction(disc,values)
        
        % Returns a matrix (or set of matrices) using the designated
        % discretization parameters.
        out = matrix(disc,varargin)
        
        % Returns a linear system RHS using the designated discretization
        % parameters.
        b = rhs(disc,f,varargin)
        
        % Reduces (projects) block rows to make space for the constraints.
        [PA,P] = reduce(disc,blocks)
        
    end
    
end
