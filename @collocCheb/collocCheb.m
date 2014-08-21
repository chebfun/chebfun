classdef collocCheb < colloc
%COLLOCCHEB   Abstract class for collocation discretization of operators.
%   COLLOCCHEB is a partial implementation of COLLOC using 1st or 2nd kind
%   Chebyshev points.
%
% See also COLLOC, COLLOC1, COLLOC2, CHEBDISCRETIZATION.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = collocCheb(varargin)
            disc = disc@colloc(varargin{:});
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS (IMPLEMENTED BY THIS ABSTRACT CLASS.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Barycentric differentiation matrix:
        D = baryDiffMat(x, w, k, t);
        
        % Dimension reduction for operator matrix.
        [PA, P, PS] = reduce(disc, A, S)
        
    end
    
end