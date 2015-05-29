classdef spherefun < separableApprox
    
    % TODO: Improve documentation of input options.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function f = spherefun(varargin)
            % The main spherefun constructor!
            
            % Return an empty CHEBFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
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
        %f = conj(f);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Hidden = true )
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PUBLIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        % The main bulk of the SPHEREFUN constructor:
        g = constructor(g, op, dom, varargin);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        %         % Convert Chebyshev coefficients to values:
        %         X = coeffs2vals(U);
        
        %         % Convert values to Chebyshev coefficients:
        %         X = vals2coeffs(U);
        
        %         % Tensor product of Chebyshev points:
        %         [xx, yy] = chebpts2(nx, ny, domain, kind);
        
        %         % Outer-product of two chebfuns:
        %         F = outerProduct(f, g);
        
        % Converts a function in spherical coordinates to one in Cartesian
        % coordinates on the sphere.
        fdf = sphf2cartf(f, lam, th, coord);
        
        % Degree l Order m spherical harmonic.
        Y = sphharm(l,m,coord)          
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private Static methods implemented by SPHEREFUN class.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = private, Static = true )
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % DOMAIN: default for the is [-pi,pi] x [0,pi].
        % which corresponds to using colatitude for the elevation angle
        % (second input argument).  Doubled-up sphere will have a domain of
        % [-pi,pi] x [-pi,pi].        
        blockDiag       % Pivot matrices used during GE
        idxPlus
        idxMinus
        pivotIndices
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Private constant properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Constant )
        alpha = 50;  % Growth factor control.
    end
    
    %
    %     methods
    %         function g = spherefun( varargin )
    %             if( nargin == 0 )
    %
    %             else
    %                 g = constructor(g , varargin{:} );  % pass to constructor.
    %             end
    %         end
    %     end
    
end