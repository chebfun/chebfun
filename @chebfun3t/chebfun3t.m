classdef chebfun3t
%CHEBFUN3T   Class for representing functions on [a,b]x[c,d]x[e,g] using 
%   tensor product approach.
%
%   This is experimental code, not for ordinary use.
% 
%   Class for approximating smooth functions defined on finite cuboids 
%   using the vanilla flavoured full tensor product approach. The output 
%   contains a tensor of coefficients of the Chebyshev expansion of the 
%   input function.
%
%   CHEBFUN3T(F) constructs a CHEBFUN3T object representing the function F on
%   [-1,1]x[-1 1]x[-1 1]. F should be a function handle, e.g.,
%   @(x,y,z) x.*y + cos(x).*z. F should in be "vectorized" in 
%   the sense that it may be evaluated at a tensor of points and returns a 
%   tensor output. 
%
%   CHEBFUN3T(F, 'eps', ep) specifies chebfun3eps to be ep.
%
%   CHEBFUN3T(F, [A B C D E G]) specifies a cube [A B]x[C D]x[E G] where the 
%   function is defined. A, B, C, D, E and G must all be finite.
%
%   Examples: 
%   f = @(x,y,z) 1./(1 + x.^2 + y.^2 + z.^2);
%   f3t = chebfun3t(f);
%   f3t = chebfun3t(f, [-2 2 1 2 1 3]);
%   f3t = chebfun3t(f,'eps', 1e-10);
%   f3t = chebfun3t(f,[-2 2 1 2 1 3],'eps', 1e-10);
%
% See also CHEBFUN3 and CHEBFUN3V.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties        
        % COEFFS: tensor of coefficients of the Chebyshev expansion of f.
		coeffs
        
        % VSCALE: max abs of tensor of samples.
		vscale
                
        % DOMAIN: cuboid of CHEBFUN3T, default is [-1,1] x [-1,1] x [-1,1].
        domain
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
            function f = chebfun3t(varargin)
            % The main CHEBFUN3T constructor!
            
            % Return an empty CHEBFUN:
            if ( (nargin == 0) || isempty(varargin{1}) )
                return
            end
            % Call the constructor. All the work is done here:
            f = constructor(f, varargin{:});
            end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = true)
        % Unfolding of a discrete tensor.
        varargout = unfold(varargin)
        
        % Reshape a discrete matrix back to a tensor.
        varargout = fold(varargin)
        
        % Tensor-matrix contractions.
        varargout = txm(varargin)        
        
        % Convert 3D values to 3D coefficients.
        coeffs3D = vals2coeffs(vals3D);
    
        % Convert 3D coefficients to 3D values.
        vals3D = coeffs2vals(coeffs3D);
    
    end
    
    methods (Access = public)
        % Retrieve and modify preferences for this class.
        varargout = subsref(f, index);
        
        % Display a CHEBFUN3T.
        varargout = display(varargin);
        
        % sampleTest for a CHEBFUN3T object.
        varargout = sampleTest(f, varargin);
        
        % Evaluate a CHEBFUN3T.
        y = feval(f, x, y, z)
        
        % Get properties of a CHEBFUN3T object.
        out = get(f, idx);
        
        % Length of a CHEBFUN3T.
        varargout = length(f);
        
        % Number of degrees of freedom needed to represent a CHEBFUN3T.
        out = ndf(f);
        
        % Definite integral of a CHEBFUN3T over its domain.
        out = sum3(f);
        
        % Norm of a CHEBFUN3T.
        out = norm(f);

        % Plotcoeffs of a CHEBFUN3T.
        out = plotcoeffs(f)

        % Addition of two CHEBFUN3T objects.
        out = plus(f, g);
        
        % Negation of a CHEBFUN3T object.
        out = uminus(f);        
        
        % Subtraction of two CHEBFUN3T objects.
        out = minus(f, g);
        
        % Pointwise multiplication of two CHEBFUN3T objects.
        out = times(f, g);
        
        % Scalar multiplication for CHEBFUN3T objects.
        out = mtimes(f, g);
        
        % Pointwise right divide of CHEBFUN3T objects.
        out = rdivide(f, g);
        
        % Pointwise left divide of CHEBFUN3T objects.
        out = ldivide(f, g);
        
        % Right divide for CHEBFUN3T objects.
        out = mrdivide(f, g);        
        
        % Pointwise power of a CHEBFUN3T.
        out = power(varargin);
        
        % Sine of a CHEBFUN3T.
        out = sin(f);

        % Cosine of a CHEBFUN3T.
        out = cos(f);
       
        % Tangent of a CHEBFUN3T.
        out = tan(f);
       
        % Tangent of a CHEBFUN3T (in degrees).
        out = tand(f);
        
        % Hyperbolic sine of a CHEBFUN3T.
        out = sinh(f);
       
        % Hyperbolic cosine of a CHEBFUN3T.
        out = cosh(f);

        % Hyperbolic tangent of a CHEBFUN3T.
        out = tanh(f);
       
        % Exponential of a CHEBFUN3T.
        out = exp(f);
              
        % Test whether a CHEBFUN3T object is empty.
        out = isempty(f);
       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Hidden = true, Static = false )
        
        % Test if two CHEBFUN3T objects have the same domain.
        out = domainCheck(f, g);
        
        % Vertical scale of a CHEBFUN3T
        out = vertscale(f);
        
    end    
    
end