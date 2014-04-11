classdef colloc < chebDiscretization
    %COLLOC   Abstract class for collocation discretization of operators.
    %
    %   See COLLOC1, COLLOC2, CHEBDISCRETIZATION.
    
    %   COLLOC is a partial implementation of CHEBDISCRETIZATION that creates
    %   scaffolding common to first-kind and second-kind points. COLLOC cannot
    %   be used directly as a discretization for linops. Both COLLOC1 and
    %   COLLOC2 are full implementations.
    
    %  Copyright 2013 by The University of Oxford and The Chebfun Developers.
    %  See http://www.chebfun.org for Chebfun information.
    
    properties (Access=private)
        % Stores LU factors of a matrix, for repeated solves at fixed size:
        mldivideData = [];
    end
    
    properties
        % Store the size of the input space relative to disc.dimension
        inputSpace = [];
    end
    
    methods
        function disc = colloc(source, dimension, domain)
            %COLLOC    Collocation discretization constructor.
            % (Called by subclasses for parts in common.)
            
            % If SOURCE is not passed, return an empty object.
            if ( isempty(source) )
                return
            end
            
            % Attach SOURCE and the DOMAIN information to the object.
            disc.source = source;
            disc.domain = source.domain;
            
            % Assign DIMENSIONS and DOMAIN if they were passed.
            if ( nargin > 1 )
                disc.dimension = dimension;
                if ( nargin > 2 )
                    disc.domain = domain;
                end
            end

            disc.inputDimension = disc.getInputDimension(source);
            
        end
        
    end
    
    % These must be implemented by a subclass.
    methods ( Abstract )
        C = cumsum(disc)    % indefinite integration
        D = diff(disc,m)    % differentiation
        % Points where function values are represented.
        [x, w] = functionPoints(disc)
        % Points where equations are enforced.
        [x, w] = equationPoints(disc)
    end
    
    methods ( Static )
        
        [x, w, v] = points(varargin);
        
        space = getInputSpace(L);
        
    end
    
end