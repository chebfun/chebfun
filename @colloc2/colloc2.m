classdef colloc2 < colloc
%COLLOC2    Collocation discretization on 2nd kind points.
%   COLLOC2 is an implementation of CHEBDISCRETIZATION that implements
%   spectral collocation on 2nd-kind Chebyshev points for differential and
%   integral operators. 
%
%   Linear algebra operations generally take O(N^3) flops, where N is
%   determined automatically to resolve the solution. You can control the
%   allowed values of N through CHEBOPPREF.
%
%   See also CHEBDISCRETIZATION, CHEBOPPREF, CHEBOP. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

    % No subclass-specific properties needed, and no special constructor
    % either. 
    
    methods
        function disc = colloc2(varargin)
            disc = disc@colloc(varargin{:});
       end
    end
    
end