classdef diskfun < separableApprox
% DISKFUN class for representing functions on the unit disk.
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
% If F is a function handle then it should allow for vectorized 
% evaluations.
%
% If F is a matrix, F = (f_ij), the numbers fij are used as function values
% at the tensor Fourier-Chebyshev points on the rectangle [-pi,pi]x[0,1].
%
% The DISKFUN software system is based on: 
%
% A. Townsend, H. Wilber, and G. Wright, Computing with functions in
% spherical and polar geometries II: The disk, submitted, 2016. 
%
% See also CHEBFUN2, SPHEREFUN, DISKFUNV.
    
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
        
        % Fast Poisson solver: 
        u = poisson( f, bc, m, n );
        
        % Converts a function in polar coordinates to one in Cartesian
        % coordinates on the disk
        fdf = pol2cartf(f, r, th);
        
        % Disk harmonics
        Y = harmonic(l,m,type) ;   
        
        %convert coeffs to a diskfun
        f = coeffs2diskfun(CFS); 
        
        
        %isCartesian: determines whether input is in
        %polar or cartesian coords. 
        iscart = isCartesian(varargin);
        
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
        %alpha = 50;  % Growth factor control.
    end
    
   
end