classdef colloc1 < colloc
%COLLOC1    Collocation discretization on 1st kind points.
%   COLLOC1 is an implementation of COLLOC that implements spectral
%   collocation on 1st-kind Chebyshev points for differential and integral
%   operators.
%
%   Linear algebra operations generally take O(N^3) flops, where N is determined
%   automatically to resolve the solution. You can control the allowed values of
%   N through CHEBOPPREF.
%
% See also COLLOC, CHEBDISCRETIZATION, CHEBOPPREF, CHEBOP. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    % No subclass-specific properties needed, and no special constructor either.

    methods
        function disc = colloc1(varargin)            
            disc = disc@colloc(varargin{:});
        end
    end
    
end
