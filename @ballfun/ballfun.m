classdef ballfun
%BALLFUN class for representing functions on the unit ball.
% 
%   Class for approximating functions defined on the unit ball. The 
%   functions should be smooth. Ballfun objects can be constructed in a 
%   variety of ways: 
%
%    1. ballfun(F), where F is a function handle in Cartesian coordinates
%       (x,y,z), e.g., @(x,y,z) x.*y.*z + cos(x).
%    2. ballfun(F, 'spherical'), where F is a function handle in spherical
%       coordinates (r,lambda,theta). Here, r is the radial variable
%       satisfying 0 <= r <= 1, lambda is the azimuthal variable and
%       satisfies -pi <= lambda <= pi and theta is the polar angle and
%       satisfies 0 <= theta < pi, e.g., 
%       @(r,lambda,theta) cos(r.*cos(lambda).*sin(theta))
%    3. ballfun(F), where F = (f_ijk) is a tensor of numbers, used as
%       function values at tensor of Chebyshev-Fourier-Fourier points in
%       the intrinsic spherical coordinate system, i.e.,
%       [0,1]x[-pi,pi]x[0,pi].
%    4. ballfun(F, 'coeffs'), where F = (f_ijk) is a tensor of Chebyshev-
%       Fourier-Fourier coefficients.
%    5. ballfun(F, [m,n,p]) returns a representation of F using a degree m
%       Chebyshev approximation of F in the r (radial) direction, a degree n
%       trigonometric approximation in the lambda (azimuthal) direction and
%       a degree p trigonometric approximation in the theta (polar) direction.
%
% If F is a function handle then it should allow for vectorized evaluations.
% 
% The BALLFUN software system is based on: 
%
% N. Boulle, and A. Townsend, Computing with functions on the ball, in
% preparation.
%
% See also CHEBFUN2, CHEBFUN3, DISKFUN, SPHEREFUN, BALLFUNV

% Copyright 2019 by The University of Oxford and The CHEBFUN Developers.
% See http://www.chebfun.org/ for CHEBFUN information.
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = ballfun(varargin)
            % The main BALLFUN constructor!
            
            % Return an empty BALLFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                f.coeffs = [];
                f.isReal = [];
                f.domain = []; 
            else
                % Call the constructor, all the work is done here:
                f = constructor(f, varargin{:});
            end
       
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % Chebyshev-Fourier-Fourier coefficients array of a BALLFUN.
        coeffs   
        
        % Boolean value designating whether the BALLFUN represents a
        % real-valued function. This allows us to always return a real result
        % for things like evaluating a BALLFUN.
        isReal
        
        %domain for non-doubled ball
        domain
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        % Convert to Chebyshev--Fourier--Fourier values
        VALS = coeffs2vals(CFS);
        
        % Convert to Chebyshev--Fourier--Fourier values
        CFS = vals2coeffs(VALS);
        
        % Compute the solid harmonics
        f = solharm(l,m,varargin);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )
        
    end
    
    methods ( Access = private, Static = false )
        
    end
end
