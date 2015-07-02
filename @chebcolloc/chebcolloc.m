classdef chebcolloc < valsDiscretization
%CHEBCOLLOC   Abstract class for collocation discretization of operators.
%   CHEBCOLLOC is a partial implementation of VALSDISCRETIZATION using 1st or
%   2nd kind Chebyshev points.
%
% See also VALSDISCRETIZATION, CHEBCOLLOC1, CHEBCOLLOC2, OPDISCRETIZATION.

%  Copyright 2015 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = chebcolloc(varargin)
            disc = disc@valsDiscretization(varargin{:});
        end
        
        % Dimension reduction for operator matrix.
        [PA, P, PS] = reduce(disc, A, S)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS (IMPLEMENTED BY THIS ABSTRACT CLASS.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Barycentric differentiation matrix:
        D = baryDiffMat(x, w, k, t);
        
    end
    
end
