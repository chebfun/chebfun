classdef (InferiorClasses = {?chebfun}) adchebfun 
%ADCHEBFUN   
%
% See also 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADCHEBFUN Class Description:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Properties of ADCHEBFUN objects.
    properties ( Access = public )
        % Domain of the ADchebfun:
        domain
        % The CHEBFUN / constant:
        func
        % The jacobian:
        jacobian
%         % Linearity information:
%         isConstant = 1;        
%         % Zero information:
%         isZero = 0;
        
    end
    
    %% CLASS CONSTRUCTOR:
    
    methods

        function obj = adchebfun(varargin)
            if ( nargin == 1 && isa(varargin{1}, 'chebfun') )
                obj.func = varargin{1};
            else
                obj.func = chebfun(varargin{:});
            end
            dom = obj.func.domain;
            obj.domain = dom;
            obj.jacobian = linop.eye(dom);
        end
        
        function f = diff(f, k)
            if ( nargin < 2 )
                k = 1; 
            end
            f.func = diff(f.func, k);
            f.jacobian = linop.diff(f.domain, k)*f.jacobian;
        end
        
        function f = times(f, g)
            % TODO: This assumes both are adchebfuns
            if ( ~isa(f, 'adchebfun') )
                g.jacobian = linop.diag(f)*g.jacobian;
                g.func = f.*g.func;
                f = g;
            elseif  ( ~isa(g, 'adchebfun') )
                f.jacobian = linop.diag(g)*f.jacobian;
                f.func = f.func.*g;
            else
                f.jacobian = linop.diag(f.func)*g.jacobian + linop.diag(g.func)*f.jacobian;
                f.func = times(f.func, g.func);
            end
        end
        
        function f = mtimes(f, g)
            if ( isnumeric(g) )
                f.func = f.func*g;
                f.jacobian = f.jacobian*g;
            elseif ( isnumeric(f) )
                g.func = f*g.func;
                g.jacobian = f*g.jacobian;
                f = g;
            else
                error
            end
        end        
        
        function f = plus(f,g)
            if ( ~isa(f, 'adchebfun') )
                f = plus(g, f);
                return
            end
            if ( ~isa(g, 'adchebfun') )
                g = adchebfun(g);
            end          
            f.func = f.func + g.func;
            f.jacobian = f.jacobian + g.jacobian;
        end
        
        function f = minus(f, g)
            f = plus(f, -g);
        end
        
        function f = uminus(f)
            f.func = -f.func;
            f.jacobian = -f.jacobian;
        end
        
        function f = sin(f)
            f.jacobian = linop.diag(cos(f.func))*f.jacobian;
            f.func = sin(f.func);
        end
        
        function f = cos(f)
            f.jacobian = linop.diag(-sin(f.func))*f.jacobian;
            f.func = cos(f.func);
        end
        
        function f = exp(f)
            f.func = exp(f.func);
            f.jacobian = linop.diag(f.func)*f.jacobian;
        end
        
        function f = power(f, b)
            % TODO: Assumes f is a CHEBFUN and b is numeric for now.
            f.jacobian = linop.diag(b*power(f.func, b-1)) * f.jacobian;
            f.func = power(f.func, b);
        end
        
        function f = sum(f)
            f.func = sum(f.func);
            f.jacobian = linop.sum(f.domain)*f.jacobian;
        end
        
        function f = feval(f, x)
            E = linop.eval(f.domain);
            f.jacobian = E(x)*f.jacobian;
            f.func = feval(f.func, x);
        end
        
        function out = get(f, prop)
            % Allow access to any of F's properties via GET.
            out = f.(prop);
        end
            
        
    end
   
    %% STATIC METHODS IMPLEMENTED BY BNDFUN CLASS.
    methods ( Static = true ) 
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
    end
end


   
