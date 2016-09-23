classdef trigspec < coeffsDiscretization
%TRIGSPEC    Fourier spectral method in coefficient space.
%   TRIGSPEC is an implementation of OPDISCRETIZATION that implements a
%   Fourier spectral method in coefficient space.
%
% See also TRIGCOLLOC.
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        coeffs        % Coefficients of the operator.
        outputSpace   % The range of the operator.
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = trigspec(varargin)
            disc = disc@coeffsDiscretization(varargin{:});
            % No dimension adjustment are required for TRIGSPEC.
            disc.dimAdjust = 0;
            disc.projOrder = 0;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        function tech = returnTech()
            %RETURNTECH    Return the appropriate tech to use for TRIGSPEC.
            tech = @trigtech;
        end
        
        % Differentiation matrices for TRIGSPEC.
        D = diffmat(N, m, flag)
        
        % Multiplication matrices for TRIGSPEC.
        D = multmat(N, f)
        
    end
    
end
