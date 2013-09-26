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
        isConstant = 1;        
        
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
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        function f = abs(f) %#ok<MANU>
            error('CHEBFUN:AD:abs:NotDifferentiable', ...
                'ABS() is not Frechet differentiable.');
        end
        
        function f = airy(k, f)
            if ( nargin == 1 )
                f = k;
                k = 0;
            end
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.diag(airy(k+1, f.func))*f.jacobian;
            f.func = airy(k, f.func);
        end
        
        function f = cos(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.diag(-sin(f.func))*f.jacobian;
            f.func = cos(f.func);
        end
        
        function f = cumsum(f)
            f.func = cumsum(f);
            f.jacobian = linop.cumsum(f.domain)*f.jacobian; 
%             f.isConstant = f.isConstant;
        end
        
        function f = diff(f, k)
            if ( nargin < 2 )
                k = 1; 
            end
            f.func = diff(f.func, k);
            f.jacobian = linop.diff(f.domain, k)*f.jacobian;          
%             f.isConstant = f.isConstant;
        end

        function f = exp(f)
            f.isConstant = iszero(f.jacobian);
            f.func = exp(f.func);
            f.jacobian = linop.diag(f.func)*f.jacobian;
%             f.domain = f.domain;            
        end
        
        function f = feval(f, x)
            E = linop.evalAt(x, f.domain);
            f.jacobian = E*f.jacobian;
            f.func = feval(f.func, x);
%             f.isConstant = f.isConstant;            
        end 
        
        function out = get(f, prop)
            % Allow access to any of F's properties via GET.
            out = f.(prop);
        end  
        
        function u = jacreset(u)
            u.jacobian = linop.eye(u.domain);
            u.isConstant = 1;
        end 
        
        function f = log(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.diag(1./f.func)*f.jacobian;
            f.func = log(f.func);
            f = updateDomain(f);
        end
        
        function f = minus(f, g)
            f = plus(f, -g);
        end
                
        function f = mtimes(f, g)
            if ( ~isa(g, 'adchebfun') )
                f.func = f.func*g;
                f.jacobian = f.jacobian*g;
            elseif ( ~isa(f, 'adchebfun') )
                g.func = f*g.func;
                g.jacobian = f*g.jacobian;
                f = g;
            else
                error('CHEBFUN:AD:mtimes:dims', ...
                    ['Matrix dimensions must agree. Use f.*g to multiply ' ...
                     'two chebfun objects.']);
            end
        end        
        
        function f = plus(f, g)
            if ( ~isa(f, 'adchebfun') )
                g.func = f + g.func;
                f = g;
            elseif ( ~isa(g, 'adchebfun') )
                f.func = f.func + g;
            else        
                f.isConstant = f.isConstant & g.isConstant;
                f.func = f.func + g.func;
                f.jacobian = f.jacobian + g.jacobian;
            end
            f = updateDomain(f);

        end
        
        function f = power(f, b)
            if ( isa(f, 'adchebfun') && isa(b, 'adchebfun') )
                f.isConstant = iszero(f.jacobian) & iszero(g.jacobian) & ...
                    f.isConstant & b.isConstant; 
                tmp = power(f.func, b.func);
                f.jacobian = diag(b.func.*f.func^(b.func-1))*f.jacobian + ...
                    diag(tmp.*log(f.func))*b.jacobian; 
                f.func = tmp;             
            elseif ( isa(f, 'adchebfun') )
                if ( isnumeric(b) )
                    if ( b == 1 )
                        % Nothing to do!
                        return
                    elseif ( b == 0 )
                        f.func = power(f.func, 0);
                        f.jacobian = 0*f.jacobian;
                        f.isConstant = 1;
                        return
                    end
                end
                f.isConstant = iszero(f.jacobian);
                f.jacobian = linop.diag(b*power(f.func, b-1))*f.jacobian;
                f.func = power(f.func, b);
                
            elseif ( isa(b, 'adchebfun') )
                b.isConstant = iszero(b.jacobian);
                b.func = power(f, b.func);
                b.jacobian = linop.diag(b.func.*log(f))*b.jacobian; 
                b.domain = b.func.domain;
                f = b;
            end
            f = updateDomain(f);
        end
      
        function u = seed(u, k, m)
            dom = u.domain;
            I = linop.eye(dom);
            Z = linop.zeros(dom);
            blocks = cell(1, m);
            for j = 1:m
                blocks{1,j} = Z;
            end
            blocks{1,k} = I;
            u.jacobian = chebmatrix(blocks);
        end        
   
        function f = sin(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.diag(cos(f.func))*f.jacobian;
            f.func = sin(f.func);
        end
        
        function f = sum(f)
            f.func = sum(f.func);
            f.jacobian = linop.sum(f.domain)*f.jacobian;
        end
        
        function f = times(f, g)
            if ( ~isa(f, 'adchebfun') )
                g.jacobian = linop.diag(f)*g.jacobian;
                g.func = f.*g.func;
%                 g.isConstant = g.isConstant;
                g = updateDomain(g);
                f = g;
            elseif  ( ~isa(g, 'adchebfun') )
                f.jacobian = linop.diag(g)*f.jacobian;
                f.func = f.func.*g;
%                 f.isConstant = f.isConstant;
            else
                f.isConstant = ...
                    ( f.isConstant & g.isConstant) & ...
                      ( ( all(iszero(f.jacobian)) || all(iszero(g.jacobian)) ) | ...
                        ( iszero(f.jacobian) & iszero(g.jacobian) ) );
                f.jacobian = linop.diag(f.func)*g.jacobian + linop.diag(g.func)*f.jacobian;
                f.func = times(f.func, g.func);
            end
            f = updateDomain(f);
        end

        function f = uminus(f)
            f.func = -f.func;
            f.jacobian = -f.jacobian;
%             f.isConstant = f.isConstant;            
        end
        
        function f = updateDomain(f)
            if ( isa(f.func, 'chebfun') )
                f.domain = union(f.domain, f.func.domain);
            end
            if ( isa(f.jacobian, 'chebmatrix') )
                f.domain = union(f.domain, domain(f.jacobian));
            else
                f.domain = union(f.domain, f.jacobian.domain);
            end
        end
        
    end
   
    %% STATIC METHODS IMPLEMENTED BY BNDFUN CLASS.
    methods ( Static = true ) 
        
    end

end


   
