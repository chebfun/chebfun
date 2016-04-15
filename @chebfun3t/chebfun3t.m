classdef chebfun3t
%CHEBFUN3T   CHEBFUN3T class for representing functions on [a,b]x[c,d]x[e,g].
% 
%   Class for approximating smooth functions defined on finite boxes using 
%   the vanilla flavoured full tensor product approach. The output contains 
%   a tensor of coefficients of the Chebyshev expansion of the input function.
%
%   CHEBFUN3T(F) constructs a CHEBFUN3T object representing the function F on
%   [-1,1]x[-1 1]x[-1 1]. F should be a function handle, e.g.,
%   @(x,y,z) x.*y + cos(x).*z. F should in be "vectorized" in 
%   the sense that it may be evaluated at a tensor of points and returns a 
%   tensor output. 
%
%   CHEBFUN3T(F, 'eps', ep) specifies chebfun3eps to be epsVal.
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties        
        % COEFFS: tensor of coefficients of the Chebyshev expansion of f
		coeffs
        
        % VSCALE: max abs of tensor of samples
		vscale
                
        % DOMAIN: box of CHEBFUN3T, default is [-1,1] x [-1,1] x [-1,1]
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
            % Call the constructor, all the work is done here:
            f = constructor(f, varargin{:});
            end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = true)
        % Tensor unfolding
        varargout = unfold(varargin)
        
        % Reshape the matrix back to a tensor
        varargout = fold(varargin)
        
        % Tensor contractions
        varargout = txm(varargin)        
        
        % vals2coeffs
        coeffs3D = vals2coeffs(vals3D);
    
        % coeffs2vals
        vals3D = coeffs2vals(coeffs3D);
    
    end
    
    methods (Access = public)
        % Retrieve and modify preferences for this class.
        varargout = subsref(f, index);
        
        % Display a CHEBFUN3T.
        varargout = display(varargin);
        
        % sampleTest
        varargout = sampleTest(f, varargin);
        
        % Evaluate a CHEBFUN3T.
        y = feval(f, x, y, z)
        
        % Get properties of a CHEBFUN3T object.
        out = get(f, idx);
        
        % Length of a CHEBFUN3tT
        varargout = length(f);
        
        % Number of degrees of freedom needed to represent a CHEBFUN3T
        out = ndf(f);

        % Definite integral of a CHEBFUN3T over its domain
        out = sum3(f);
        
        out = plotcoeffs(f)

        % Addition of two CHEBFUN3T objects
        out = plus(f, g);
        
        % Negation of a CHEBFUN3T object
        out = uminus(f);        
        
        % Subtraction of two CHEBFUN3T objects
        out = minus(f, g);
        
        % Multiplication of two CHEBFUN3T objects
        %out = times(f, g);
        out = times(f, g, tol);
        
        out = power(varargin);
        
       out = sin(f);

       %Cosine of a CHEBFUN3T
       out = cos(f);
       
       out = tan(f);
       
       out = tand(f);
       
       out = tanh(f);
       
       out = exp(f);
       
       out = sinh(f);
       
       out = cosh(f);
    end
    
end
