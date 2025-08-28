classdef diskfun < separableApprox
%DISKFUN   DISKFUN class for representing functions on the unit disk.
% 
%   Class for approximating smooth functions defined on the unit disk. The 
%   functions should be smooth. Diskfun objects can be constructed in a 
%   variety of ways: 
%
%    1. Y = diskfun(F), where F is a function handle in Cartesian 
%       coordinates (x,y), e.g., @(x,y) x.*y + cos(x).
%    2. Y = diskfun(F, 'polar'), where F is a function handle in polar 
%       coordinates (theta,r). Here, theta is the angular variable satisfying
%       -pi <= theta <= pi, and r is the radial variable satisfying 0 <= r < 1,
%       e.g., @(theta,r) cos(r.*sin(theta))
%    3. Y = diskfun(F), where F is a matrix of numbers.
%
% If F is a function handle then it should allow for vectorized evaluations.
%
% If F is a matrix, F = (f_ij), the numbers fij are used as function values
% at the tensor Fourier-Chebyshev points on the rectangle [-pi,pi]x[0,1].
%
% DISKFUN(F, k) returns a rank k approximation to F.
%
% DISKFUN(F, [m n]) returns a representation of F using a degree m
% Chebyshev approximation of F in the radial direction and a degree n
% trigonometric approximation in the angular direction. The result is
% compressed in low rank form and the rank k is still determined adaptively
% (satisfying k<=min(m,n)+1).
% 
% The DISKFUN software system is based on: 
%
% H. Wilber, A. Townsend, and G. Wright, Computing with functions in
% spherical and polar geometries II: The disk, SIAM. J. Sci. Comput., 39-4
% (2017), C238-C262.
%
% See also CHEBFUN2, SPHEREFUN, DISKFUNV.

% Copyright 2017 by The University of Oxford and The CHEBFUN Developers.
% See http://www.chebfun.org/ for CHEBFUN information.

% TODO: Include documentation of fixed eps construction

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = diskfun(varargin)
            % The main diskfun constructor!
            
            % Return an empty DISKFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                f.domain = [-pi pi 0 1];
                return
            end
            % Call the constructor, all the work is done here:
            f = constructor(f, varargin{:});       
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
  
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Hidden = true )
        % Make f have BMC-II symmetry.
        f = projectOntoBMCII(f);
        
    end
      
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PUBLIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
     
        % The main bulk of the DISKFUN constructor:
        g = constructor(g, op, dom, varargin);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Fast Poisson & Helmholtz solvers: 
        u = poisson( f, bc, m, n );
        u = helmholtz( f, k, bc, m, n); 
        
        % Converts a function in polar coordinates to one in Cartesian
        % coordinates on the disk
        fdf = pol2cartf(f, r, th);
        
        % Disk harmonics
        Y = harmonic(l,m,type);
        
        % Convert coeffs to a diskfun
        f = coeffs2diskfun(CFS);    
        
        
        varargout = coeffs2vals(U, varargin); 
        varargout = vals2coeffs(U, varargin);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private Static methods implemented by DISKFUN class.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )
    

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = public)
        % DOMAIN: default is [-pi,pi] x [0,1] which corresponds to using 
        % polar coordinates. 
        % Doubled-up disk will have a domain of [-pi,pi] x [-1,1].
        idxPlus
        idxMinus
        nonZeroPoles = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private constant properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Constant )
        % TODO: Add support for this constant here.
        %alpha = 50;  % Growth factor control.
    end
    
   
end