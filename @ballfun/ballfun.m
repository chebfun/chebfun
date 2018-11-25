classdef ballfun
%BALLFUN class for representing functions on the unit ball.
% 
%   Class for approximating functions defined on the unit ball. The 
%   functions should be smooth.
%
% BALLFUN(F) constructs a BALLFUN object representing the function F on
% the unit ball. F can have the following form:
%    1. A function handle in (x,y,z), e.g., @(x,y,z) x.*y.*z + cos(x).
%    2. A function handle in spherical coordinates (r,lambda,theta), where
%       lambda is the azimuthal variable and satisfies -pi <= lambda <= pi
%       and theta is the polar angle and satisfies 0 <= theta < pi,
%       e.g., @(r,lambda,theta) cos(r.*cos(lambda).*sin(theta))
%    3. A matrix of numbers. 
% If F is a function handle then it should allow for vectorized evaluations.
%
% If F is a matrix, F = (f_ijk), the numbers fijk are used as function values
% at tensor of Chebyshev-Fourier-Fourier points in the intrinsic spherical 
% coordinate system, i.e., [0,1]x[-pi,pi]x[0,pi].
% 
% The BALLFUN software system is based on: 
%
% N. Boulle, and A. Townsend, Computing with functions in the ball, in
% preparation.
%
% See also CHEBFUN2, DISKFUN, SPHEREFUN, BALLFUNV

% Copyright 2018 by The University of Oxford and The CHEBFUN Developers.
% See http://www.chebfun.org/ for CHEBFUN information.
    
    properties
        
        % Chebyshev-Fourier-Fourier coefficients array of a BALLFUN.
        coeffs   
        
        % Boolean value designating whether the BALLFUN represents a
        % real-valued function. This allows us to always return a real result
        % for things like evaluating a BALLFUN.
        isReal   
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = ballfun(varargin)
            % The main BALLFUN constructor!
            
            % Return an empty BALLFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            
            % Call the constructor, all the work is done here:
            f = constructor(f, varargin{:});
       
        end
    end
    
    methods ( Access = public, Static = true )
        
        % Convert to Chebyshev--Fourier--Fourier values
        VALS = coeffs2vals(CFS);
        
        % Convert to Chebyshev--Fourier--Fourier values
        CFS = vals2coeffs(VALS);
        
        % Compute the solid harmonics
        f = solharm(l,m);
                
    end
    
    methods ( Access = private, Static = true )
        
    end
    
    methods ( Access = private, Static = false )
        
    end
    
end
