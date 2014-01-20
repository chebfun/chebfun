classdef (InferiorClasses = {?chebfun}) adchebfun 
%ADCHEBFUN   A class consisting of a CHEBFUN and derivative information.
%
% See also CHEBFUN, LINOP, CHEBOP.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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
        % Linearity information:
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
            obj.jacobian = linBlock.eye(dom);
        end
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        function f = abs(f) %#ok<MANU>
            error('CHEBFUN:AD:abs:NotDifferentiable', ...
                'ABS() is not Frechet differentiable.');
        end
        
        function f = airy(k, f)
            % F = AIRY(K,F)   Airy function of an ADCHEBFUN.
            %
            % F = AIRY(F) where F is a CHEBFUN is the same as above, with K = 0,
            
            % Default value of K
            if ( nargin == 1 )
                f = k;
                k = 0;
            end
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(airy(k + 1, f.func))*f.jacobian;
            % Update CHEBFUN part.
            f.func = airy(k, f.func);
        end

        function f = acos(f)
            % F = ACOS(F)   ACOS of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-1./sqrt(1 - f.func.^2))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acos(f.func);
        end
        
        function f = acosd(f)
            % F = ACOSD(F)   ACOSD of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-(180/pi)./sqrt(1 - f.func.^2))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acosd(f.func);
        end
        
        function f = acosh(f)
            % F = ACOSH(F)   ACOSH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./sqrt(f.func.^2 - 1))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acosh(f.func);
        end
        
        function f = acot(f)
            % F = ACOT(F)   ACOT of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-1./(1 + f.func.^2))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acot(f.func);
        end
        
        function f = acotd(f)
            % F = ACOTD(F)   ACOTD of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-(180/pi)./(1 + f.func.^2))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acotd(f.func);
        end
        
        function f = acoth(f)
            % F = ACOTH(F)   ACOTH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-1./(f.func.^2 - 1))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acoth(f.func);
        end
        
        function f = acsc(f)
            % F = ACSC(F)   ACSC of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-1./(abs(f.func).*sqrt(f.func.^2 - 1)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acsc(f.func);
        end
        
        function f = acscd(f)
            % F = ACSCD(F)   ACSCD of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-(180/pi)./(abs(f.func).*sqrt(f.func.^2 - 1)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acscd(f.func);
        end
                
        function f = acsch(f)
            % F = ACSCH(F)   ACSCH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-1./(f.func.*sqrt(1 + f.func.^2)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acsch(f.func);
        end
        
        function f = asec(f)
            % F = ASEC(F)   ASEC of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./(abs(f.func).*sqrt(f.func.^2-1)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = asec(f.func);
        end
                        
        function f = asecd(f)
            % F = ASECD(F)   ASECD of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult((180/pi)./(abs(f.func).*sqrt(f.func.^2-1)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = asecd(f.func);
        end
                        
        function f = asech(f)
            % F = ASECH(F)   ASECH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-1./(f.func.*sqrt(1-f.func.^2)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = asech(f.func);
        end
                        
        function f = asin(f)
            % F = ASIN(F)   ASIN of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./sqrt(1-f.func.^2))*f.jacobian;
            % Update CHEBFUN part.
            f.func = asin(f.func);
        end
                        
        function f = asind(f)
            % F = ASIND(F)   ASIND of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult((180/pi)./sqrt(1-f.func.^2))*f.jacobian;
            % Update CHEBFUN part.
            f.func = asind(f.func);
        end
                        
        function f = asinh(f)
            % F = ASINH(F)   ASINH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./sqrt(f.func.^2+1))*f.jacobian;
            % Update CHEBFUN part.
            f.func = asinh(f.func);
        end
                                
        function f = atan(f)
            % F = ATAN(F)   ATAN of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./(1+f.func.^2))*f.jacobian;
            % Update CHEBFUN part.
            f.func = atan(f.func);
        end
                                
        function f = atand(f)
            % F = ATAND(F)   ATAND of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult((180/pi)./(1+f.func.^2))*f.jacobian;
            % Update CHEBFUN part.
            f.func = atand(f.func);
        end
                                
        function f = atanh(f)
            % F = ATANH(F)   ATANH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./(1-f.func.^2))*f.jacobian;
            % Update CHEBFUN part.
            f.func = atanh(f.func);
        end
        
        function g = besselj(nu, f)
            % G = BESSELJ(NU, F)    Bessel-J function of an ADCHEBFUN.
            %
            % See also chebfun/besselj.
            
            % Initialise an empty ADCHEBFUN
            g = adchebfun;
            % Copy the domain information
            g.domain = f.domain;
            % Linearity information
            g.isConstant = iszero(f.jacobian);
            % Function composition
            g.func = besselj(nu, f.func);
            % Derivative computation
            g.jacobian = operatorBlock.mult(-besselj(nu+1,f.func) + ...
                nu*(g.func)./f.func)*f.jacobian;
        end
        
        function f = cos(f)
            % F = COS(F)   COS of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linBlock.mult(-sin(f.func))*f.jacobian;
            % Update CHEBFUN part.
            f.func = cos(f.func);
        end
        
        function f = cosd(f)
            % F = COSD(F)   COSD of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-pi/180*sind(f.func))*f.jacobian;
            % Update CHEBFUN part.
            f.func = cosd(f.func);
        end
        
        function f = cosh(f)
            % F = COSH(F)   COSH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(sinh(f.func))*f.jacobian;
            % Update CHEBFUN part.
            f.func = cosh(f.func);
        end
                
        function f = cot(f)
            % F = COT(F)   COT of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-csc(f.func).^2)*f.jacobian;
            % Update CHEBFUN part.
            f.func = cot(f.func);
        end
                
        function f = cotd(f)
            % F = COTD(F)   COTD of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-(pi/180)*cscd(f.func).^2)*f.jacobian;
            % Update CHEBFUN part.
            f.func = cotd(f.func);
        end
        
        function f = coth(f)
            % F = COTH(F)   COTH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-csch(f.func).^2)*f.jacobian;
            % Update CHEBFUN part.
            f.func = coth(f.func);
        end
        
        function f = cumsum(f, k)
            % F = CUMSUM(F, K)   CUMSUM of an ADCHEBFUN
            
            % By default, compute first anti-derivative
            if ( nargin < 2 )
                k = 1;
            end
            
            % Update CHEBFUN part
            f.func = cumsum(f, k);
            % Update derivative part
            f.jacobian = linop.cumsum(f.domain, k)*f.jacobian;
            % CUMSUM is a linear operation, so no need to update linearity info.
            % f.isConstant = f.isConstant;
        end
        
        function f = diff(f, k)
            % F = DIFF(F, K)   DIFF of an ADCHEBFUN
            
            % By default, compute first derivative
            if ( nargin < 2 )
                k = 1; 
            end
            
            % Update CHEBFUN part
            f.func = diff(f.func, k);
            % Update derivative part
            f.jacobian = linop.diff(f.domain, k)*f.jacobian;
            % DIFF is a linear operation, so no need to update linearity info.
            % f.isConstant = f.isConstant;
        end
       
        function f = erf(f)
            % F = ERF(F)   ERF of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(2*exp(-f.func.^2)/sqrt(pi))*f.jacobian;
            % Update CHEBFUN part
            f.func = erf(f.func);      
        end
                
        function f = erfc(f)
            % F = ERFC(F)   ERFC of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-2*exp(-f.func.^2)/sqrt(pi))*f.jacobian;
            % Update CHEBFUN part
            f.func = erfc(f.func);      
        end
        
        function f = erfcinv(f)
            % F = ERFCINV(F)   ERFCINV of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part
            f.func = erfcinv(f.func);
            % Update derivative part
            f.jacobian = operatorBlock.mult(exp(f.func.^2)*sqrt(pi)/2)*f.jacobian;        
        end
        
        function g = erfcx(f)
            % F = ERFCX(F)   ERFCX of an ADCHEBFUN.
            
            % Need to copy F to G, as we need info about both functions to
            % compute derivatives below.
            g = f;
            % Linearity information
            g.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part
            g.func = erfcx(f.func);
            % Update derivative part
            g.jacobian = operatorBlock.mult(-2/sqrt(pi) + 2*(f.func).*(g.func))*f.jacobian;        
        end
        
        function f = erfinv(f)
            % F = ERFINV(F)   ERFINV of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part
            f.func = erfinv(f.func);
            % Update derivative part
            f.jacobian = operatorBlock.mult(exp(f.func.^2)*sqrt(pi)/2)*f.jacobian;
        end
        
        function f = exp(f)
            % F = EXP(F)   EXP of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part
            f.func = exp(f.func);
            % Update derivative part
            f.jacobian = linBlock.mult(f.func)*f.jacobian;        
        end
        
        function f = feval(f, x)
            % F = FEVAL(F,X)    Evaluate an ADCHEBFUN F at point X.
            %
            % The output will be an ADCHEBFUN with derivative representing
            % evaluation at the point X (and the derivative of the input, as
            % dictated by the chain rule).
            
            % Create an feval linear operator at the point X.
            E = linop.feval(x, f.domain);
            % Update derivative part
            f.jacobian = E*f.jacobian;
            % Update CHEBFUN part
            f.func = feval(f.func, x);
            % Evaluation is a linear operation, so no need to update linearity
            % information.
            % f.isConstant = f.isConstant;
        end 
        
        function out = get(f, prop, pos)
            % Allow access to any of F's properties via GET.
            if nargin == 2
                out = vertcat(f.(prop));
            else
                out = f(pos).(prop);
            end
        end
        
        function out = getElement(f, pos)
            % Return the pos-th element
            out = f(pos);
        end
        
        function u = jacreset(u)
            u.jacobian = linop.eye(u.domain);
            u.isConstant = 1;
        end 
        
        function f = log(f)
            % F = LOG(F)   LOG of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult(1./f.func)*f.jacobian;
            f.func = log(f.func);
            f = updateDomain(f);
        end
                
        function f = log1p(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult(1./(f.func + 1))*f.jacobian;
            f.func = log1p(f.func);
            f = updateDomain(f);
        end
 
        function f = log2(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult(1./(log(2)*f.func))*f.jacobian;
            f.func = log2(f.func);
            f = updateDomain(f);
        end
        
        function f = log10(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult(1./(log(10)*f.func))*f.jacobian;
            f.func = log10(f.func);
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
                f.jacobian = operatorBlock.mult(b*power(f.func, b-1))*f.jacobian;
                f.func = power(f.func, b);
                
            elseif ( isa(b, 'adchebfun') )
                b.isConstant = iszero(b.jacobian);
                b.func = power(f, b.func);
                b.jacobian = operatorBlock.mult(b.func.*log(f))*b.jacobian; 
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
   
        function g = sec(f)
            g = f;
            g.isConstant = iszero(f.jacobian);
            g.func = sec(f.func);
            g.jacobian = operatorBlock.mult(tan(f.func).*g.func)*f.jacobian;
        end
           
        function g = secd(f)
            g = f;
            g.isConstant = iszero(f.jacobian);
            g.func = secd(f.func);
            g.jacobian = operatorBlock.mult(pi/180*tand(f.func).*g.func)*f.jacobian;
        end
        
        function g = sech(f)
            g = f;
            g.isConstant = iszero(f.jacobian);
            g.func = sech(f.func);
            g.jacobian = operatorBlock.mult(-tanh(f.func).*g.func)*f.jacobian;
        end
        
        function f = sin(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = linBlock.mult(cos(f.func))*f.jacobian;
            f.func = sin(f.func);
        end
        
        function f = sinh(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult(cosh(f.func))*f.jacobian;
            f.func = sinh(f.func);
        end

        function out = subsref(f, index)
            switch index(1).type
                case '()'
                    out = feval(f, index.subs{1});
                case '.'
%                     out = vertcat(f.(index(1).subs));
                    out = f.(index(1).subs);
            end
        end
        
        function f = sum(f)
            f.func = sum(f.func);
            f.jacobian = linop.sum(f.domain)*f.jacobian;
        end
        
        function f = tan(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult(sec(f.func).^2)*f.jacobian;
            f.func = tan(f.func);
        end
        
        function f = tand(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult((pi/180)*secd(f.func).^2)*f.jacobian;
            f.func = tand(f.func);
        end
        
        function f = tanh(f)
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult(sech(f.func).^2)*f.jacobian;
            f.func = tanh(f.func);
        end
        
        function f = times(f, g)
            if ( isnumeric(f) ) || ( isnumeric(g) )
                f = mtimes(f, g);
            elseif ( ~isa(f, 'adchebfun') )
                g.jacobian = linBlock.mult(f)*g.jacobian;
                g.func = f.*g.func;
%                 g.isConstant = g.isConstant;
                g = updateDomain(g);
                f = g;
            elseif  ( ~isa(g, 'adchebfun') )
                f.jacobian = linBlock.mult(g)*f.jacobian;
                f.func = f.func.*g;
                f = updateDomain(f);
%                 f.isConstant = f.isConstant;
            else
                f.isConstant = ...
                    ( f.isConstant & g.isConstant) & ...
                      ( ( all(iszero(f.jacobian)) || all(iszero(g.jacobian)) ) | ...
                        ( iszero(f.jacobian) & iszero(g.jacobian) ) );
                f.jacobian = linBlock.mult(f.func)*g.jacobian + linBlock.mult(g.func)*f.jacobian;
                f.func = times(f.func, g.func);
                f = updateDomain(f);
            end

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


   
