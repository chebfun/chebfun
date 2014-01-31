classdef colloc2 < chebDiscretization
%COLLOC2   Collocation based on 2nd-kind Chebyshev points.
%   COLLOC2 is an implementation of CHEBDISCRETIZATION that applies spectral
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