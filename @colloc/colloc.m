classdef colloc < chebDiscretization
%COLLOC   Abstract class for collocation discretization of operators.
%   COLLOC is a partial implementation of CHEBDISCRETIZATION that creates
%   scaffolding common to first-kind and second-kind points. COLLOC cannot
%   be used directly as a discretization for linops. Both COLLOC1 and
%   COLLOC2 are full implementations.
%
%   See also COLLOC1, COLLOC2, CHEBDISCRETIZATION.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.
    
    properties ( Access = private )
        % Stores LU factors of a matrix, for repeated solves at fixed size:
        mldivideData = [];
    end
    
    methods
        function disc = colloc(source, dimension, domain)
            %COLLOC    Collocation discretization constructor.
            %   (Called by subclasses for parts in common.)
            
            if ( (nargin == 0) || isempty(source) )
                % Construct an empty COLLOC.
                return
            end
            
            % Attach SOURCE and the DOMAIN information to the object:
            disc.source = source;
            disc.domain = source.domain;
            % Determine the dimension adjustment:
            disc.dimAdjust = colloc.getDimAdjust(source);
            
            % Assign DIMENSIONS and DOMAIN if they were passed:
            if ( nargin > 1 )
                disc.dimension = dimension;
            end
            if ( nargin > 2 )
                disc.domain = domain;
            end
                
        end
        
    end
    
    % These must be implemented by a subclass.
    methods ( Abstract )
        
        % Indefinite integration:
        C = cumsum(disc)    
        
        % Differentiation:
        D = diff(disc, m)  
        
        % Points where function values are represented:
        [x, w] = functionPoints(disc)
        
        % Points where equations are enforced:
        [x, w] = equationPoints(disc)
        
    end
    
    methods ( Static )
        
        % Discretization points: (used by both colloc1 and colloc2)
        [x, w, v] = points(varargin);

    end
    
end
