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
            f.jacobian = linop.mult(airy(k + 1, f.func))*f.jacobian;
            f.func = airy(k, f.func);
        end

        function f = acos(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-1./sqrt(1 - f.func.^2))*f.jacobian;
            f.func = acos(f.func);
        end
        
        function f = acosd(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-(180/pi)./sqrt(1 - f.func.^2))*f.jacobian;
            f.func = acosd(f.func);
        end
        
        function f = acosh(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(1./sqrt(f.func.^2 - 1))*f.jacobian;
            f.func = acosh(f.func);
        end
        
        function f = acot(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-1./(1 + f.func.^2))*f.jacobian;
            f.func = acot(f.func);
        end
        
        function f = acotd(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-(180/pi)./(1 + f.func.^2))*f.jacobian;
            f.func = acotd(f.func);
        end
        
        function f = acoth(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-1./(f.func.^2 - 1))*f.jacobian;
            f.func = acoth(f.func);
        end
        
        function f = acsc(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-1./(abs(f.func).*sqrt(f.func.^2 - 1)))*f.jacobian;
            f.func = acsc(f.func);
        end
        
        function f = acscd(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-(180/pi)./(abs(f.func).*sqrt(f.func.^2 - 1)))*f.jacobian;
            f.func = acscd(f.func);
        end
                
        function f = acsch(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-1./(f.func.*sqrt(1 + f.func.^2)))*f.jacobian;
            f.func = acsch(f.func);
        end
        
        function f = asec(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(1./(abs(f.func).*sqrt(f.func.^2-1)))*f.jacobian;
            f.func = asec(f.func);
        end
                        
        function f = asecd(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult((180/pi)./(abs(f.func).*sqrt(f.func.^2-1)))*f.jacobian;
            f.func = asecd(f.func);
        end
                        
        function f = asech(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-1./(f.func.*sqrt(1-f.func.^2)))*f.jacobian;
            f.func = asech(f.func);
        end
                        
        function f = asin(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(1./sqrt(1-f.func.^2))*f.jacobian;
            f.func = asin(f.func);
        end
                        
        function f = asind(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult((180/pi)./sqrt(1-f.func.^2))*f.jacobian;
            f.func = asind(f.func);
        end
                        
        function f = asinh(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(1./sqrt(f.func.^2+1))*f.jacobian;
            f.func = asinh(f.func);
        end
                                
        function f = atan(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(1./(1+f.func.^2))*f.jacobian;
            f.func = atan(f.func);
        end
                                
        function f = atand(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult((180/pi)./(1+f.func.^2))*f.jacobian;
            f.func = atand(f.func);
        end
                                
        function f = atanh(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(1./(1-f.func.^2))*f.jacobian;
            f.func = atanh(f.func);
        end
        
        function f = cos(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linop.mult(-sin(f.func))*f.jacobian;
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
            f.jacobian = linop.mult(f.func)*f.jacobian;
%             f.domain = f.domain;            
        end
        
        function f = feval(f, x)
            E = linop.feval(x, f.domain);
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
            f.jacobian = linop.mult(1./f.func)*f.jacobian;
            f.func = log(f.func);
            f = updateDomain(f);
        end
        
        function f = minus(f, g)
            f = plus(f, -g);
        end
                
        function [normF, normLoc] = norm(f, varargin)
            if nargout == 2
                [normF, normLoc]  = norm(f.func, varargin{:});
            else
                normF = norm(f.func, varargin{:});
            end
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
                f.jacobian = linop.mult(b*power(f.func, b-1))*f.jacobian;
                f.func = power(f.func, b);
                
            elseif ( isa(b, 'adchebfun') )
                b.isConstant = iszero(b.jacobian);
                b.func = power(f, b.func);
                b.jacobian = linop.mult(b.func.*log(f))*b.jacobian; 
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
            f.jacobian = linop.mult(cos(f.func))*f.jacobian;
            f.func = sin(f.func);
        end
        
        function f = sum(f)
            f.func = sum(f.func);
            f.jacobian = linop.sum(f.domain)*f.jacobian;
        end
        
        function f = times(f, g)
            if ( ~isa(f, 'adchebfun') )
                g.jacobian = linop.mult(f)*g.jacobian;
                g.func = f.*g.func;
%                 g.isConstant = g.isConstant;
                g = updateDomain(g);
                f = g;
            elseif  ( ~isa(g, 'adchebfun') )
                f.jacobian = linop.mult(g)*f.jacobian;
                f.func = f.func.*g;
%                 f.isConstant = f.isConstant;
            else
                f.isConstant = ...
                    ( f.isConstant & g.isConstant) & ...
                      ( ( all(iszero(f.jacobian)) || all(iszero(g.jacobian)) ) | ...
                        ( iszero(f.jacobian) & iszero(g.jacobian) ) );
                f.jacobian = linop.mult(f.func)*g.jacobian + linop.mult(g.func)*f.jacobian;
                f.func = times(f.func, g.func);
            end
            f = updateDomain(f);
        end

        function f = uminus(f)
            f.func = -f.func;
            f.jacobian = -f.jacobian;
%             f.isConstant = f.isConstant;            
        end
        
%         function f = vertcat(varargin)
%             if ( nargin > 1 )
%                 f = chebmatrix(varargin.');
%             else
%                 f = varargin{1};
%             end
%         end
        
        function f = updateDomain(f)
            if ( isa(f.func, 'chebfun') )
                f.domain = union(f.domain, f.func.domain);
            end
            f.domain = union(f.domain, f.jacobian.domain);
        end
        
    end
   
    %% STATIC METHODS IMPLEMENTED BY BNDFUN CLASS.
    methods ( Static = true ) 
        
    end

end


   
