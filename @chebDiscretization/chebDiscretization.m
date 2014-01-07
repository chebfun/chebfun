classdef (Abstract) chebDiscretization 
    
    % Objects of this class store a source (either a linBlock or a linop),
    % domain, and discretization dimension. Calling the matrix() method causes
    % the source to be discretized at the relevant parameters. In the linop
    % case, this includes the imposition of any side and continuity conditions. 
    
    properties
        domain = []
        dimension = []
        source = []
    end
        
    properties (Dependent)
        numIntervals
    end

    methods
        function n = get.numIntervals(disc)
            n = length(disc.domain) - 1;
        end        

        function t = isempty(disc)
            t = isempty(disc.source);
        end
        
        function t = isFactored(disc)
            % This method gives a discretization a chance to store matrix
            % factors for the purpose of short-circuiting the linsolve process.
            % By default it never happens. 
            t = false;
        end
        
        % Each discretization is free to replace standard backslash with
        % whatever it likes for the discrete linear system solution. This is the
        % default case. 
        function [x,disc] = mldivide(disc,A,b)
            x = A\b;
        end           
        
    end
        
    methods (Abstract)
        values = toValues(disc,f)
        f = toFunction(disc,values)
        A = matrix(disc,varargin)
        b = rhs(disc,f,varargin)
    end
    
end
