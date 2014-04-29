classdef colloc2 < chebDiscretization
%COLLOC2   Collocation based on 2nd-kind Chebyshev points.
%   COLLOC2 is an implementation of CHEBDISCRETIZATION that applies spectral
%   collocation using 2nd kind Chebyshev points for differential and integral
%   operators and systems. 
%
%   The default discretization type to use is set by CHEBOPPREF. You can also
%   use CHEBOPPREF to create a preferences object and change its
%   'discretization' property. 
%
%   Linear algebra operations with COLLOC2 operators generally take O(N^3)
%   flops, where in most contexts N is determined automatically to resolve the
%   solution. The allowed values of N are governed by the 'dimensionValues'
%   property in CHEBOPPREF. You can also set the maximum N (including systems
%   and piecewise definitions) through the 'maxTotalLength' property.
%
%   See also CHEBOPPREF, CHEBDISCRETIZATION.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

    properties (Access=private)
        mldivideData = [];  % stores LU factors of a matrix for repeated solves
    end
     
    methods
        function disc = colloc2(source, dimension, domain)
        % COLLOC2 constructor
        
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
        
end
