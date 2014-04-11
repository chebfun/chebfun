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
    
    methods
        function disc = colloc(source, dimension, domain)
            %COLLOC    Collocation discretization constructor.
            
            % Called by subclasses for parts in common. 
            % If SOURCE is not passed, return an empty object.
            if isempty(source)
                return
            end
            % Attach SOURCE and the DOMAIN information to the object.
            disc.source = source;
            disc.domain = source.domain;
            
            % Assign DIMENSIONS and DOMAIN if they were passed.
            if ( nargin > 1 )
                disc.dimension = dimension;
                if nargin > 2
                    disc.domain = domain;
                end
            end
        end
        
    end
        
    % These must be implemented by a subclass.
    methods ( Abstract )
        C = cumsum(disc)    % indefinite integration
        D = diff(disc,m)    % differentiation 
        % Points where function values are represented.
        [x,w] = functionPoints(disc)
        % Points where equations are enforced. 
        [x,w] = equationPoints(disc)
    end
    
    methods ( Static )
        function [x, w, v] = points(varargin)
%POINTS    Discretization points.
%   X = COLLOC.POINTS(DISC,KIND) returns KIND-kind points using the domain and
%   dimension stored in DISC. KIND must be either 1 or 2.
%
%   An alternate calling sequence is COLLOC.POINTS(DOMAIN,DIMENSION,KIND).
%
%   [X, W, V] = COLLOC.POINTS(DISC,KIND) also returns Clenshaw-Curtis
%   weights and barycentric interpolation weights, respectively.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

            if ( nargin == 2 )
                disc = varargin{1};
                kind = varargin{2};
                d = disc.domain;
                n = disc.dimension;

            elseif ( nargin == 3 )
                d = varargin{1};
                n = varargin{2};
                kind = varargin{3};
                
            else
                error('Must be called with two or three arguments.')
                
            end
            
            numInt = length(n);
            
            % Create output as cells for convenience.
            x = cell(numInt, 1);
            w = cell(1, numInt);
            v = cell(numInt, 1);
            for k = 1:numInt
                
                % To save time, don't call again unless the dimension has changed.
                if ( k==1 ) || ( n(k) ~= n(k-1) )
                    if ( kind == 2 )
                        [x0, w0, v0] = chebtech2.chebpts(n(k));
                    else
                        [x0, w0, v0] = chebtech1.chebpts(n(k));
                    end
                end
                
                % The points and weights returned by the CHEBTECH methods above live on
                % [-1, 1]. Transform the points to the interval we're working on.
                dif = (d(k + 1) - d(k))/2;
                x{k} = x0*dif + (d(k +1) + d(k))/2;
                
                if ( nargout > 1 )
                    w{k} = w0*dif;
                    if ( nargout > 2 )
                        v{k} = v0;
                    end
                end
                
            end
            
            % Convert output to vectors.
            x = cell2mat(x);
            if ( nargout > 1 )
                w = cell2mat(w);
                if ( nargout > 2 )
                    v = cell2mat(v);
                end
            end
            
        end
        
    end
    
end