classdef (InferiorClasses = {?chebfun}) adchebfun 
%ADCHEBFUN   A class for supporting automatic differentiation in Chebfun.
%
%   The ADCHEBFUN class allows Chebfun to compute Frechet derivatives of
%   nonlinear operators. It also supports linearity detection.
% 
%   This class is not intended to be called directly by the end user.
%
% See also CHEBFUN, LINBLOCK, LINOP, CHEBOP.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADCHEBFUN class description:
%
% The ADCHEBFUN class is used by the CHEBOP class to supply Frechet derivatives
% for solving nonlinear boundary-value problems (BVPs) of ordinary differential
% equations (ODEs). It also enables the system to determine whether operators
% are linear or not, which allows a convenient for syntax for specifying linear
% operators, and accessing methods specific to linear operators, such as EIGS()
% and EXPM().
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Properties of ADCHEBFUN objects.
    properties ( Access = public )
        % FUNC: A CHEBFUN (or double) , which corresponds to the function the
        % ADCHEBFUN  represents.
        func
        
        % JACOBIAN: The Frechet derivative of the function the ADCHEBFUN
        % represents, with respect to a selected basis variable. The basis
        % variable is selected at the start of computation with ADCHEBFUN
        % objects, and has the identity operators as its Frechet derivative.
        jacobian
        
        % ISCONSTANT: This is used for linearity detection. A value equal to 1
        % indicates that the function the ADCHEBFUN represents has Frechet
        % derivatives which has a constant value with respect to the selected
        % basis variable.
        isConstant = 1;  
        
        % DOMAIN: Domain of the ADchebfun.
        domain
    end
    
    %% CLASS CONSTRUCTOR:
    
    methods
        
        function obj = adchebfun(varargin)
            % TODO: Document. Do we want to allow passing derivatives here for
            % instant seeding at construction time?
            if ( nargin == 1 && isa(varargin{1}, 'chebfun') )
                obj.func = varargin{1};
            else
                obj.func = chebfun(varargin{:});
            end
            dom = obj.func.domain;
            obj.domain = dom;
            obj.jacobian = operatorBlock.eye(dom);
        end
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        function f = abs(f) %#ok<MANU>
            % ABS   ABS is not Frechet differentiable, so an error is thrown.
            error('CHEBFUN:AD:abs:NotDifferentiable', ...
                'ABS() is not Frechet differentiable.');
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
            f.jacobian = operatorBlock.mult(-(180/pi)./sqrt(1 - f.func.^2))*...
                f.jacobian;
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
            f.jacobian = operatorBlock.mult(-(180/pi)./(1 + f.func.^2))*...
                f.jacobian;
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
            f.jacobian = operatorBlock.mult(-1./(abs(f.func).* ...
                sqrt(f.func.^2 - 1)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acsc(f.func);
        end
        
        function f = acscd(f)
            % F = ACSCD(F)   ACSCD of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-(180/pi)./(abs(f.func).* ...
                sqrt(f.func.^2 - 1)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acscd(f.func);
        end
                
        function f = acsch(f)
            % F = ACSCH(F)   ACSCH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-1./(f.func.* ...
                sqrt(1 + f.func.^2)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = acsch(f.func);
        end
        
                
        function f = airy(k, f)
            % F = AIRY(K,F)   Airy function of an ADCHEBFUN.
            
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
        
        function f = asec(f)
            % F = ASEC(F)   ASEC of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./(abs(f.func).* ...
                sqrt(f.func.^2-1)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = asec(f.func);
        end
                        
        function f = asecd(f)
            % F = ASECD(F)   ASECD of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult((180/pi)./(abs(f.func).* ...
                sqrt(f.func.^2-1)))*f.jacobian;
            % Update CHEBFUN part.
            f.func = asecd(f.func);
        end
                        
        function f = asech(f)
            % F = ASECH(F)   ASECH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-1./(f.func.* ...
                sqrt(1-f.func.^2)))*f.jacobian;
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
            f.jacobian = operatorBlock.mult((180/pi)./sqrt(1-f.func.^2))* ...
                f.jacobian;
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
        
        function f = besselj(nu, f)
            % BESSELJ(NU, F)   Bessel-J function of an ADCHEBFUN.

            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Function composition
            tmp = besselj(nu, f.func);
            % Derivative computation
            f.jacobian = operatorBlock.mult(-besselj(nu+1, f.func) + ...
                nu*tmp./f.func)*f.jacobian;
            f.func = tmp;
        end
        
        function f = cos(f)
            % F = COS(F)   COS of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult(-sin(f.func))*f.jacobian;
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
            f.jacobian = operatorBlock.mult(-(pi/180)*cscd(f.func).^2)* ...
                f.jacobian;
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
        
        function g = csc(f)
            % F = CSC(F)   CSC of an ADCHEBFUN.
            
            % Copy F to the output G to enable reuse of computed function value
            g = f;
            % Linearity information
            g.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part.
            g.func = csc(f.func);
            % Update derivative part
            g.jacobian = operatorBlock.mult(-cot(f.func).*g.func)*f.jacobian;
        end
        
        function g = cscd(f)
            % F = CSCD(F)   CSCD of an ADCHEBFUN.
            
            % Copy F to the output G to enable reuse of computed function value
            g = f;
            % Linearity information
            g.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part.
            g.func = cscd(f.func);
            % Update derivative part
            g.jacobian = operatorBlock.mult(-pi/180*cotd(f.func).*g.func)*f.jacobian;
        end  
        
        function g = csch(f)
            % F = CSCH(F)   CSCH of an ADCHEBFUN.
            
            % Copy F to the output G to enable reuse of computed function value
            g = f;
            % Linearity information
            g.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part.
            g.func = csch(f.func);
            % Update derivative part
            g.jacobian = operatorBlock.mult(-coth(f.func).*g.func)*f.jacobian; 
        end       
        
        function f = cumprod(f)
            % F = CUMPROD(F)    CUMPROD of an ADCHEBFUN
            %
            % See also CHEBFUN/CUMPROD().
            
            f = exp(cumsum(log(f)));
        end
        
        function f = cumsum(f, k)
            % F = CUMSUM(F, K)   CUMSUM of an ADCHEBFUN
            
            % By default, compute first anti-derivative
            if ( nargin < 2 )
                k = 1;
            end
            
            % Update CHEBFUN part
            f.func = cumsum(f.func, k);
            % Update derivative part
            f.jacobian = operatorBlock.cumsum(f.domain, k)*f.jacobian;
            % CUMSUM is a linear operation, so no need to update linearity info.
            % f.isConstant = f.isConstant;
        end
        
        function f = diff(f, k)
            % F = DIFF(F, K)   DIFF of an ADCHEBFUN
            
            % By default, compute first derivative:
            if ( nargin < 2 )
                k = 1; 
            end
            
            % Update CHEBFUN part:
            f.func = diff(f.func, k);
            % Update derivative part:
            f.jacobian = operatorBlock.diff(f.domain, k)*f.jacobian;
            % DIFF is a linear operation, so no need to update linearity info.
            % f.isConstant = f.isConstant;
        end
        
        function [sm, cm, dm] = ellipj(f, m)
            % [SM, CM, DM] = ELLIPJ(F, M)    Ellip-J function of an ADCHEBFUN.
            
            % Either F or M could be an ADCHEBFUN, but be we currently only
            % support the case where the first argument is an ADCHEBFUN
            if ( ~isa(f, 'adchebfun') )
               error('CHEBFUN:ADCHEBFUN:ellipj', ...
                   ['Currently, ADCHEBFUN only supports calls to ELLIPJ() ' ...
                   'where the first input is a ADCHEBFUN.']);
            end
            
            % Copy F to the output G to enable reuse of computed function value:
            sm = f;
            % Linearity information:
            sm.isConstant = iszero(f.jacobian);
            
            % Copy the ADCHEBFUN SM to CM and DM also:
            cm = sm;
            dm = sm;
            
            % Do the function computation. Again, either F or M could have been
            % the ADCHEBFUN:
            [smtemp, cmtemp, dmtemp] = ellipj(f.func, m);
            
            % Assign the functional part to SM:
            sm.func = smtemp;
            
            % We know we always want the derivative info about SM:
            sm.jacobian = operatorBlock.mult(cmtemp.*dmtemp)*f.jacobian;
            
            % Compute as much derivative information is required, depending on
            % the number of outputs requested:
            if ( nargout >= 2)
                
                % Assign the function part to CM:
                cm.func = cmtemp;
                % Derivative computation:
                cm.jacobian = operatorBlock.mult(-smtemp.*dmtemp);
                
                if ( nargout > 2 )  % Also want DM!
                    
                    % Function part:
                    dm.func = dmtemp;
                    % Derivative:
                    dm.jacobian = operatorBlock.mult(-m.*smtemp.*cmtemp);
                    
                end
                
            end
            
        end
       
        function f = erf(f)
            % F = ERF(F)   ERF of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(2*exp(-f.func.^2)/sqrt(pi))* ...
                f.jacobian;
            % Update CHEBFUN part
            f.func = erf(f.func);      
        end
                
        function f = erfc(f)
            % F = ERFC(F)   ERFC of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(-2*exp(-f.func.^2)/sqrt(pi))* ...
                f.jacobian;
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
            f.jacobian = operatorBlock.mult(-exp(f.func.^2)*sqrt(pi)/2)* ...
                f.jacobian;        
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
            g.jacobian = operatorBlock.mult(-2/sqrt(pi) + ...
                2*(f.func).*(g.func))*f.jacobian;        
        end
        
        function f = erfinv(f)
            % F = ERFINV(F)   ERFINV of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part
            f.func = erfinv(f.func);
            % Update derivative part
            f.jacobian = operatorBlock.mult(exp(f.func.^2)*sqrt(pi)/2)* ...
                f.jacobian;
        end
        
        function f = exp(f)
            % F = EXP(F)   EXP of an ADCHEBFUN.

            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part
            f.func = exp(f.func);
            % Update derivative part
            f.jacobian = operatorBlock.mult(f.func)*f.jacobian;        
        end
        
        function f = expm1(f)
            % F = EXPM1(F)   EXPM1 of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(exp(f.func))*f.jacobian;  
            % Update CHEBFUN part
            f.func = expm1(f.func);      
        end
        
        function f = feval(f, x)
            % F = FEVAL(F,X)    Evaluate an ADCHEBFUN F at point X.
            %
            % The output will be an ADCHEBFUN with derivative representing
            % evaluation at the point X (and the derivative of the input, as
            % dictated by the chain rule).
            
            % Create an feval linear operator at the point X.
            E = functionalBlock.feval(x, f.domain);
            % Update derivative part
            f.jacobian = E*f.jacobian;
            % Update CHEBFUN part
            f.func = feval(f.func, x);
            % Evaluation is a linear operation, so no need to update linearity
            % information.
            % f.isConstant = f.isConstant;
        end 
        
        function f = fred(K, f, varargin)
            % FRED   Fredholm operator.
            
            % Update CHEBFUN part
            f.func = fred(K, f.func);
            
            % Update derivative part
            f.jacobian = operatorBlock.fred(f.domain, K, varargin{:})*f.jacobian;
            
            % FRED is a linear operation, so no need to update linearity info.
            % f.isConstant = f.isConstant;
            
        end
        
        function out = get(f, prop, pos)
            %GET   Access properties of ADCHEBFUN objects.
            %   P = GET(F, PROP) returns the property P specified in the string
            %   PROP from the CHEBFUN F. Valid entries for the string PROP are:
            %       'DOMAIN'     -   The domain of definintion of F.
            %       'FUNC'       -   The chebfun/scalar part of F.
            %       'JACOBIAN'   -   The derivative of F w.r.t. the seeding
            %                        variable.
            %       'ISCONSTANT' -   Whether F has a constant derivative w.r.t.
            %                        the seeding variable (useful for linearity
            %                        detection).
            %
            %   In case F is an ADCHEBFUN array, use the call
            %       P = GET(F, PROP, POS)
            %   to obtain the desired property of the element in the position
            %   POS only.
           
            % Allow access to any of F's properties via GET.
            if nargin == 2
                out = vertcat(f.(prop));
            else
                out = f(pos).(prop);
            end
        end
        
        function out = getElement(f, pos)
            % TODO: Document
            
            % Return the pos-th element
            out = f(pos);
        end
        
        function f = heaviside(f) %#ok<MANU>
            % HEAVISIDE is not Frechet differentiable, so an error is thrown.
            error('CHEBFUN:AD:heaviside:NotDifferentiable', ...
                'HEAVISIDE() is not Frechet differentiable.');
        end
        
        function out = hscale(f)
            % HSCALE    Horizontal scale of the FUNC part of an ADCHEBFUN.
            %
            % See also: CHEBFUN/HSCALE
            out = hscale(f.func);
        end
        
        function varargout = integral(varargin)
            % F = INTEGRAL(F)   INTEGRAL(F) synonym for SUM(F).
            
            [varargout{1:nargout}] = sum(varargin{:});
        end
        
        function u = jacreset(u)
            % U = JACRESET(U)
            % TODO: Document
            % TODO: Dimension?
            u.jacobian = operatorBlock.eye(u.domain);
            u.isConstant = 1;
        end 
        
        function u = jump(u, x, c)
            % U = JUMP(U)       JUMP of an ADCHEBFUN
            %
            % 
            if ( nargin < 3)
                c = 0;
            end
            
            % Update the domain, introducing a break at the jump location
            u.domain = union(u.domain, x);
            % Compute the value of the jump
            u.func = jump(u.func, x, c);
            % Derivative part
            u.jacobian = functionalBlock.jump(x, u.domain, 0)*u.jacobian;
        end 
        
        function l = length(f)
            % LENGTH(F) where F is an ADCHEBFUN is the same as LENGTH(F.FUNC)
            l = length(f.func);
        end
        
        function f = log(f)
            % F = LOG(F)   LOG of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./f.func)*f.jacobian;
            % Update CHEBFUN part
            f.func = log(f.func);
            % Need to update domain in case breakpoints were introduced
            f = updateDomain(f);
        end
                
        function f = log1p(f)
            % F = LOG1P(F)   LOG1P of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./(f.func + 1))*f.jacobian;
            % Update CHEBFUN part
            f.func = log1p(f.func);
            % Need to update domain in case breakpoints were introduced
            f = updateDomain(f);
        end
 
        function f = log2(f)
            % F = LOG2(F)   LOG2 of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./(log(2)*f.func))*f.jacobian;
            % Update CHEBFUN part
            f.func = log2(f.func);
            % Need to update domain in case breakpoints were introduced
            f = updateDomain(f);
        end
        
        function f = log10(f)
            % F = LOG10(F)   LOG10 of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(1./(log(10)*f.func))*f.jacobian;
            % Update CHEBFUN part
            f.func = log10(f.func);
            % Need to update domain in case breakpoints were introduced
            f = updateDomain(f);
        end
        
        function varargout = loglog(f, varargin)
            % LOGLOG    log-log plot of the CHEBFUN part of an ADCHEBFUN
            [varargout{1:nargout}] = loglog(f.func, varargin{:});
        end
        
        function f = minus(f, g)
            % -     Subtraction of ADCHEBFUN objects
            f = plus(f, -g);
        end
        
        function f = mrdivide(f, g)
            %/    Right matrix divide for CHEBFUN objects.
            %     F/A divides the ADCHEBFUN F by the scalar A.
            
            % If second input is numeric, call rdivide. Otherwise, throw an error.
            if ( isnumeric(g) )
                f = rdivide(f, g);
            else
                error('CHEBFUN:AD:mrdivide:dims', ...
                    ['Matrix dimensions must agree. Use f./g to divide' ...
                    'ADCHEBFUN objects.']);
            end
        end
        
        function f = mtimes(f, g)
            %*   ADCHEBFUN multiplication.
            %     A*F and F*A multiplies the ADCHEBFUN F by the scalar A.
   
            % If either input is numeric, call times. Otherwise, throw an error.
            if ( isnumeric(f) || isnumeric(g) )
                f = times(f, g);  
            else
                error('CHEBFUN:AD:mtimes:dims', ...
                    ['Matrix dimensions must agree. Use f.*g to multiply ' ...
                     'two ADCHEBFUN or CHEBFUN objects.']);
            end
        end   
                
        function varargout = norm(f, varargin)
            % NORM(F, K)    Norm of ADCHEBFUN objects.
            %
            % Input argument follow the expected pattern from CHEBFUN/norm.
            %
            % See also CHEBFUN/norm.
            
            % TODO: Do we want this method to return an ADCHEBFUN? Makes sense
            % in the 2-norm case, in particular for 2-norm squared.
            [varargout{1:nargout}] = norm(f.func, varargin{:});
        end     
        
        function varargout = plot(f, varargin)
            % PLOT      Plot the CHEBFUN part of an ADCHEBFUN
            [varargout{1:nargout}] = plot(f.func, varargin{:});
        end
        
        function f = plus(f, g)
            % +     Addition of ADCHEBFUN objects
            
            % If F is not an ADCHEBFUN, we know G is, so add to the CHEBFUN part
            % of G.
            if ( ~isa(f, 'adchebfun') )
                g.func = f + g.func;
                f = g;      % Swap for output argument
            
            % If G is not an ADCHEBFUN, we know F is, so add to the CHEBFUN part
            % of F.
            elseif ( ~isa(g, 'adchebfun') )
                f.func = f.func + g;
            
            % ADCHEBFUN + ADCHEBFUN
            else
                % Update linearity information
                f.isConstant = f.isConstant & g.isConstant;
                % Derivative part
                f.func = f.func + g.func;
                % Derivative part
                f.jacobian = f.jacobian + g.jacobian;
            end
            
            % Need to update domain in case new breakpoints were introduced
            f = updateDomain(f);

        end
        
        
        function f = pow2(f)
            % F = POW2(F)   POW2 of an ADCHEBFUN
            
            f = power(2, f);
        end        
        
        function f = power(f, b)
            %.^   ADCHEBFUN power
            %
            % See also: chebfun/power
            
            % ADCHEBFUN.^ADCHEBFUN
            if ( isa(f, 'adchebfun') && isa(b, 'adchebfun') )
                % Linearity information
                f.isConstant = iszero(f.jacobian) & iszero(b.jacobian) & ...
                    f.isConstant & b.isConstant;
                % Temporarily store the function value to be returned                
                tmp = power(f.func, b.func);
                % Derivative information
                f.jacobian = ...
                    operatorBlock.mult(b.func.*f.func.^(b.func-1))*f.jacobian + ...
                    operatorBlock.mult(tmp.*log(f.func))*b.jacobian; 
                % Assign the function value
                f.func = tmp;             
                
            % ADCHEBFUN.^SCALAR or ADCHEBFUN.^CHEBFUN    
            elseif ( isa(f, 'adchebfun') )
                if ( isnumeric(b) )
                    if ( b == 1 )
                        % Nothing to do!
                        return
                    elseif ( b == 0 )
                        f.func = power(f.func, 0);
                        f.jacobian = 0*f.jacobian;
                        f.isConstant = true(size(f.jacobian));
                        return
                    end
                end
                % Linearity information
                f.isConstant = iszero(f.jacobian);
                % Update derivative of function value
                f.jacobian = operatorBlock.mult(b.*power(f.func, b-1))*f.jacobian;
                f.func = power(f.func, b);
                
            % SCALAR.^ADCHEBFUN or CHEBFUN.^ADCHEBFUN
            elseif ( isa(b, 'adchebfun') )
                b.isConstant = iszero(b.jacobian);
                b.func = power(f, b.func);
                b.jacobian = operatorBlock.mult(b.func.*log(f))*b.jacobian; 
                % Swap variables to get output of method
                f = b;
            end
            
            % Check whether new breakpoints have to be introduced.
            f = updateDomain(f);
        end
      
        function f = prod(f)
            % F = PROD(F)       PROD of an ADCHEBFUN
            % 
            % See also chebfun/prod.
            
            f = exp(sum(log(f)));
        end     
        
        function f = rdivide(f, g)
            % ./    ADCHEBFUN division
            %
            % F./G divides F and G, where F and G may be ADHCEBFUN or CHEBFUN
            % objects or scalars.
            
            % TODO: Can the scalar and chebfun cases be merged once we can pass
            % domains to the operatorBlock.mult method?
            % ADCHEBFUN./SCALAR or ADCHEBFUN./CHEBFUN
            if ( isnumeric(g) || ~isa(g, 'adchebfun') ) 
                f = f.*(1./g);
            elseif ( isnumeric(f) )             % SCALAR./ADCHEBFUN
                % Temporariy store the function of 1./g
                tmp = 1./g.func;
                % Update function part
                g.func = f.*tmp;
                % Linearity information
                g.isConstant = iszero(g.jacobian);
                % Derivative part
                g.jacobian = operatorBlock.mult(-f.*tmp.^2)*g.jacobian;
                % Update domain in case new breakpoints were introduced
                g = updateDomain(g);
                % Swap variables for output
                f = g;
            elseif ( ~isa(f, 'adchebfun') )     % CHEBFUN./ADCHEBFUN
                % Temporariy store the function of 1./g
                tmp = 1./g.func;
                % Update function part
                g.func = f.*tmp;
                % Linearity information
                g.isConstant = iszero(g.jacobian);
                % Derivative part
                g.jacobian = operatorBlock.mult(-f.*tmp.^2)*g.jacobian;
                % Update domain in case new breakpoints were introduced
                g = updateDomain(g);
                % Swap variables for output
                f = g;
            else                                % ADCHEBFUN.*ADCHEBFUN
                % Temporarily store the function of 1./g
                tmp = 1./g.func;
                % Derivative part
                f.jacobian = operatorBlock.mult(tmp)*f.jacobian - ...
                             operatorBlock.mult(f.func.*tmp.^2)*g.jacobian;
                % Function part
                f.func = f.func.*tmp;
                
                % Rather complicated linearity information.
                % TODO: Fix this! We were doing linearity checking incorrectly
                % for ./ in the past.
                f.isConstant = iszero(g.jacobian);
                
                % Update domain in case new breakpoints were introduced.
                f = updateDomain(f);
            end
            
        end
        
        function u = seed(u, k, m)
            % TODO: Document
            dom = u.domain;
            I = operatorBlock.eye(dom);
            Z = operatorBlock.zeros(dom);
            blocks = cell(1, m);
            for j = 1:m
                blocks{1,j} = Z;
            end
            blocks{1,k} = I;
            u.jacobian = chebmatrix(blocks);
        end        
   
        function g = sec(f)
            % F = SEC(F)   SEC of an ADCHEBFUN.
            
            % Need to copy F to G, as we need info about both functions to
            % compute derivatives below.
            g = f;
            % Linearity information
            g.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part
            g.func = sec(f.func);
            % Update derivative part
            g.jacobian = operatorBlock.mult(tan(f.func).*g.func)*f.jacobian;
        end
           
        function g = secd(f)
            % F = SECD(F)   SECD of an ADCHEBFUN.
            
            % Need to copy F to G, as we need info about both functions to
            % compute derivatives below.
            g = f;
            % Linearity information
            g.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part
            g.func = secd(f.func);
            % Update derivative part
            g.jacobian = operatorBlock.mult(pi/180*tand(f.func).*g.func)* ...
                f.jacobian;
        end
        
        function g = sech(f)
            % F = SECH(F)   SECH of an ADCHEBFUN.
            
            % Need to copy F to G, as we need info about both functions to
            % compute derivatives below.
            g = f;
            % Linearity information
            g.isConstant = iszero(f.jacobian);
            % Update CHEBFUN part
            g.func = sech(f.func);
            % Update derivative part
            g.jacobian = operatorBlock.mult(-tanh(f.func).*g.func)*f.jacobian;
        end
        
        function varargout = semilogx(f, varargin)
            % SEMILOGX   Semilogx plot of the CHEBFUN part of an ADCHEBFUN
            [varargout{1:nargout}] = semilogx(f.func, varargin{:});
        end
        
        function varargout = semilogy(f, varargin)
            % SEMILOGY   Semilogy plot of the CHEBFUN part of an ADCHEBFUN
            [varargout{1:nargout}] = semilogy(f.func, varargin{:});
        end
        
        function f = sin(f)
            % F = SIN(F)   SIN of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(cos(f.func))*f.jacobian;
            % Update CHEBFUN part
            f.func = sin(f.func);
        end
        
        function f = sinc(f)
            % F = SINC(F)  SINC of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            Jop = @(u) (pi*u.*cos(pi*u) - sin(pi*u))./(pi*u.^2);
            f.jacobian = operatorBlock.mult(compose(f.func, Jop))*f.jacobian;
            % Update CHEBFUN part
            f.func = sinc(f.func);
            
        end
        
        function f = sind(f)
            % F = SIND(F)   SIND of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(pi/180*cosd(f.func))*f.jacobian;
            % Update CHEBFUN part.
            f.func = sind(f.func);
        end
        
        function f = sinh(f)
            % F = SINH(F)   SINH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(cosh(f.func))*f.jacobian;
            % Update CHEBFUN part
            f.func = sinh(f.func);
        end
        
        function varargout = size(f, varargin)
            % SIZE(F) where F is an ADCHEBFUN is the same as SIZE(F.FUNC)
            [varargout{1:nargout}] = size(f.func, varargin{:});
        end
        
        function f = sqrt(f)
            % F = SQRT(F)   SQRT of an ADCHEBFUN
            
            f = power(f, 0.5);
        end
        
        function out = subsref(f, index)
            % TODO: Document
            switch index(1).type
                case '()'
                    out = feval(f, index.subs{1});
                case '.'
%                     out = vertcat(f.(index(1).subs));
                    out = f.(index(1).subs);
            end
        end
        
        function f = sum(f)
            % F = SUM(F)   Definite integral of an ADCHEBFUN.
            
            % TODO: SUM(F, A, B) ?
            
            % Compute the definite integral of the CHEBFUN part. This will be a
            % scalar, but we want to return and ADCHEBFUN in order to be able to
            % return the derivative information as well.
            f.func = sum(f.func);
            % Update derivative information
            f.jacobian = functionalBlock.sum(f.domain)*f.jacobian;
        end
        
        function f = tan(f)
            % F = TAN(F)   TAN of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult(sec(f.func).^2)*f.jacobian;
            % Update CHEBFUN part
            f.func = tan(f.func);
        end
        
        function f = tand(f)
            % F = TAND(F)   TAND of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            % Update derivative part
            f.jacobian = operatorBlock.mult((pi/180)*secd(f.func).^2)* ...
                f.jacobian;
            % Update CHEBFUN part
            f.func = tand(f.func);
        end
        
        function f = tanh(f)
            % F = TANH(F)   TANH of an ADCHEBFUN.
            
            % Linearity information
            f.isConstant = iszero(f.jacobian);
            f.jacobian = operatorBlock.mult(sech(f.func).^2)*f.jacobian;
            % Update CHEBFUN part
            f.func = tanh(f.func);
        end
        
        function f = times(f, g)
            % .*    ADCHEBFUN multiplication
            %
            % F.*G multiplies F and G, where F and G may be ADHCEBFUN or CHEBFUN
            % objects or scalars.
            
            % TODO: Can the scalar and chebfun cases be merged once we can pass
            % domains to the operatorBlock.mult method?
            
            if ( isnumeric(g) )                 % ADCHEBFUN.*SCALAR
                f.func = f.func*g;
                f.jacobian = f.jacobian*g;
                
            elseif ( isnumeric(f) )             % SCALAR.*ADCHEBFUN
                g.func = f*g.func;
                g.jacobian = f*g.jacobian;
                % Swap variables for output
                f = g;
            elseif  ( ~isa(g, 'adchebfun') )    % ADCHEBFUN.*CHEBFUN
                f.jacobian = operatorBlock.mult(g)*f.jacobian;
                f.func = f.func.*g;
                % Update domain in case new breakpoints were introduced.
                f = updateDomain(f);
            elseif ( ~isa(f, 'adchebfun') )     % CHEBFUN.*ADCHEBFUN
                g.jacobian = operatorBlock.mult(f)*g.jacobian;
                g.func = f.*g.func;
                % Update domain in case new breakpoints were introduced.
                g = updateDomain(g);
                % Swap variables for output
                f = g;
            else                                % ADCHEBFUN.*ADCHEBFUN
                f.jacobian = operatorBlock.mult(f.func)*g.jacobian + ...
                    operatorBlock.mult(g.func)*f.jacobian;
                f.func = times(f.func, g.func);
                
                % Rather complicated linearity information
                f.isConstant = ...
                    ( f.isConstant & g.isConstant) & ...
                    ( ( all(iszero(f.jacobian)) || all(iszero(g.jacobian)) ) | ...
                    ( iszero(f.jacobian) & iszero(g.jacobian) ) );
                
                % Update domain in case new breakpoints were introduced.
                f = updateDomain(f);
            end

        end
        
        function f = uminus(f)
            % -  Unary minus of an ADCHEBFUN
            
            % Do the obvious things...
            f.func = -f.func;
            f.jacobian = -f.jacobian;         
        end
        
        function f = uplus(f)
            % -  Unary plus of an ADCHEBFUN
            
            % This method does nothing.
        end
        
        function f = volt(K, f, varargin)
            % VOLT   Volterra operator.
            
            % Update CHEBFUN part
            f.func = volt(K, f.func);
            
            % Update derivative part
            f.jacobian = operatorBlock.volt(f.domain, K, varargin{:})*f.jacobian;
            
            % VOLT is a linear operation, so no need to update linearity info.
            % f.isConstant = f.isConstant;
            
        end
        
        function out = vscale(f)
            % VSCALE    Vertical scale of the FUNC part of an ADCHEBFUN.
            %
            % See also: CHEBFUN/VSCALE
            out = vscale(f.func);
        end
        
    end
    
    methods (Static = true)
        % Taylor testing for correctness of derivatives
        [order1, order2, nDiff2] = taylorTesting(f, hMax, numOut, plotting)
        
        % Taylor testing for correctness of derivatives of binary operators
        [order1, order2, nDiff2] = taylorTestingBinary(f, hMax, plotting)
        
        % Value testing for correctness of computed function
        error = valueTesting(f, numOut)
        
        % Value testing for correctness of computed function for binary
        % operators.
        error = valueTestingBinary(f)
    end
    
    methods ( Access = private )
        
        function f = updateDomain(f)
            % UPDATEDOMAIN      Update the domain of an ADCHEBFUN
            %
            % Various ADCHEBFUN method can cause new breakpoints to be
            % introduced in the CHEBFUN part of the ADCHEBFUN. This method
            % updates the breakpoint information in the domain field of the
            % ADCHEBFUN at the end of such methods to ensure they agree.
            
            % If the func part is a CHEBFUN
            if ( isa(f.func, 'chebfun') )
                f.domain = union(f.domain, f.func.domain);
                
                % If the func part is a scalar.
            end
            f.domain = union(f.domain, f.jacobian.domain);
        end
    end
    
end


   
