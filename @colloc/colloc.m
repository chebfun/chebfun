classdef colloc < chebDiscretization
%COLLOC   Collocation discretization of operators.
%   COLLOC is an implementation of CHEBDISCRETIZATION that applies spectral
%   collocation using 2nd kind Chebyshev points for differential and integral
%   operators and systems. To use COLLOC2 for a linop L, set
%
%     L.prefs.discretization = @colloc2;
%
%   Linear algebra operations with COLLOC2 operators generally take O(N^3)
%   flops, where N is determined automatically to resolve the solution. You can
%   control the allowed values of N through the setting
%
%     L.prefs.dimensionValues = [ values ]
%
%   If you give a single value here, the discretization will be of fixed size.
%   Note that the matrix might not use that value exactly due to breakpoints and
%   multiple variables. You can also set the maximum N through
%
%      L.prefs.maxTotalLength = N
%
%   which limits N summed over all variables; i.e. the actual matrix being
%   manipulated.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    
    properties (Access=private)
        mldivideData = [];  % stores LU factors of a matrix for repeated solves
    end
    
    methods
        function disc = colloc(source, dimension, domain)
            %COLLOC    Collocation discretization constructor.
            
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
        
    methods ( Abstract )
        C = cumsum(disc)    % indefinite integration
        D = diff(disc,m)    % differentiation 
        % Points where function values are represented.
        [x,w] = functionPoints(disc)
        % Points where equations are enforced. 
        [x,w] = equationPoints(disc)
    end
    
    methods ( Static )
        function [x, w, v] = points(disc, kind)
%POINTS    Discretization points.
%   X = COLLOC.POINTS(DISC,KIND) returns KIND-kind points using the domain and
%   dimension stored in DISC. KIND must be either 1 or 2.
%
%   [X, W] = COLLOC.POINTS(DISC,KIND) also returns Clenshaw-Curtis weights.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
            
            % Obtain useful info
            d = disc.domain;
            numInt = disc.numIntervals;
            n = disc.dimension;
            
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