classdef trigcolloc < valsDiscretization
%TRIGCOLLOC   Collocation discretization on equally spaced points.
%
%   TRIGCOLLOC is an implementation of VALSDISCRETIZATION that implements
%   spectral collocation on equi-spaced points for differential and integral
%   operators.
%
%   Linear algebra operations generally take O(N^3) flops, where N is determined
%   automatically to resolve the solution. You can control the allowed values of
%   N through CHEBOPPREF.
%
% See also COLLOC, OPDISCRETIZATION, CHEBOPPREF, CHEBOP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = trigcolloc(varargin)
            disc = disc@valsDiscretization(varargin{:});
            % No dimension adjustment are required for TRIGCOLLOC.
            disc.dimAdjust = 0; 
            disc.projOrder = 0; 
        end
        
        % Dimension reduction for operator matrix.
        [PA, P, PS] = reduce(disc, A, S);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        function tech = returnTech()
        %RETURNTECH    Return the appropriate tech to use for TRIGCOLLOC.
            tech = @trigtech;
        end
        
        % Trigonometric Fourier differentiation matrix.
        D = diffmat(N, m);
        
        % Trigonometric Fourier integration matrix.
        Q = cumsummat(N);
    
    end

end
