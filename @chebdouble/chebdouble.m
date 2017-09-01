classdef (InferiorClasses = {?chebfun}) chebdouble
%CHEBDOUBLE  Chebyshev double class. For example, DIFF means Chebyshev difference.
%
%   PDE15S likes to work with doubles (for speed). However, the problem is that
%   a call to PDE15s of the form pdeFun = @(u) diff(u) will simply call the
%   built-in DIFF method and compute a finite difference method, rather than the
%   derivative of the Chebyshev interpolant to the data.
%
%   To get around this, we use this CHEBDOUBLE class, which behaves in the same
%   way as a double for almost all methods, except for DIFF, SUM, CUMSUM, and
%   FEVAL, in which the stored values are presumbed to be values on a Chebyshev
%   grid, and the appropriate action is taken.
%
%   This class in intended solely as a worker-class for PDESOLVE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        
        % VALUES: Values of a Chebyshev interpolant.
        values = [];
        
        % DOMAIN: Domain of the interpolant.
        domain = [-1, 1];
        
        % DIFFORDER: Used to determine the highest order spatial derivative in a
        % PDE or system of PDEs.
        diffOrder = 0;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function obj = chebdouble(vals, dom)
            obj.values = vals;
            if ( nargin > 1 )
                obj.domain = dom;
            end
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DIFF  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = diff(u, k)
            %DIFF   Compute the k-th derivative of u using Chebyshev
            % differentiation matrices defined by diffmat.
            
            % Store the diffmat D as a persistent variable to allow speeding up
            % if we work with the same discretization at multiple time steps.
            % Note that the matrix D is independent of the domain, since it is
            % scaled separately below.
            persistent D
            
            % Assume first-order derivative
            if ( nargin == 1 )
                k = 1;
            end
            
            N = length(u.values);
            
            % Construct D if we don't match a previous discretization.
            if ( isempty(D) || numel(D) < k || size(D{k}, 1) ~= N )
                D{k} = chebcolloc2.diffmat(N, k); % Diffmat
            end
            
            % Interval scaling
            c = 2/diff(u.domain);     
            
            % Muliplying by the kth-order differentiation matrix
            u.values = c^k*(D{k}*u.values);
            
            % Update the difforder:
            u.diffOrder = u.diffOrder + k;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SUM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The differential operators
        function I = sum(u, a, b)
            %SIM  Compute the integral of u using Clenshaw-Curtis nodes and
            %     weights (which are stored for speed).
            
            persistent W
            
            if ( isempty(W) )
                W = {};
            end
            
            % Extract the data:
            N = length(u.values);
            
            % Deal with the 3 args case. This can be integrating a sub-domain or
            % indefinite integration. (Or integrating the whole domain...)
            if ( nargin == 3 )
                x = chebpts(N, u.domain);
                if ( length(b) > 1 )
                    if ( ~all(b == x) )
                        error('CHEBFUN:CHEBDOUBLE:sum:sumb', ...
                            ['Limits in sum must be scalars or the ', ...
                            'independent space variable (typically ''x'').']);
                    elseif ( a < x(1) )
                        error('CHEBFUN:CHEBDOUBLE:sum:sumint', ...
                            'Limits of integration outside of domain.');
                    end
                    I = cumsum(u);
                    I = I - feval(I, a);
                    return
                elseif ( length(a) > 1 )
                    if ( ~all(a == x) )
                        error('CHEBFUN:CHEBDOUBLE:sum:suma', ...
                            ['Limits in sum must be scalars or the ', ...
                            'independent space variable (typically ''x'').']);
                    elseif ( b > x(end) )
                        error('CHEBFUN:CHEBDOUBLE:sum:sumint', ...
                            'Limits of integration outside of domain.');
                    end
                    I = cumsum(u);
                    I = feval(I, b) - I;
                    return
                elseif ( a ~= x(1) || b ~= x(end) )
                    if ( a < x(1) || b > x(end) )
                        error('CHEBFUN:CHEBDOUBLE:sum:sumint', ...
                            'Limits of integration outside of domain.');
                    end
                    I = cumsum(u);
                    I = feval(I, b) - feval(I, a);
                    return
                else
                    % Sum(u, a, b) is the same as below!
                end
            end
            
            % Retrieve or compute weights::
            if ( N > 5 && numel(W) >= N && ~isempty(W{N}) )
                % Weights are already in storage!
            else
                c = diff(u.domain)/2; % Interval scaling.
                W{N} = c*chebtech2.quadwts(N);
            end
            
            % Find the sum by muliplying by the weights vector:
            I = W{N}*u;
        end
        
        function I = integral(varargin)
            I = sum(varargin{:});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CUMSUM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The differential operators
        function u = cumsum(u)
            %CUMSUM   Compute the indefinite integral of the Chebyshev
            %         interpolant to u.
            
            persistent C

            % Extract the data:
            N = length(u.values);
            c = diff(u.domain)/2; % Interval scaling.
            
            % Compute cumsum matrix:
            if ( size(C, 2) ~= N )
                C = chebcolloc2.cumsummat(N);
            end
            
            % Compute the indefinite integral:
            u.values = c*(C*u.values);
            
            % Update the difforder:
            u.diffOrder = u.diffOrder - 1;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FRED  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = fred(K, u)
            %FRED  Fredholm operator with kernel K.
            %   FRED(K, U) computes the action of the Fredholm operator with
            %   kernel K on the Chebyshev interpolant to the points in the
            %   vector U.
            
            persistent X W
            if ( isempty(W) )
                X = {};
                W = {};
            end
            
            % Extract the data:
            N = length(u);
            
            % Retrieve or compute weights::
            if ( N > 5 && numel(W) >= N && ~isempty(W{N}) )
                % Weights are already in storage!
            else
                [X{N}, W{N}] = chebpts(N, u.domain);
            end
            
            % The Fredholm operator:
            [xx, yy] = ndgrid(X{N});
            u.values = K(xx, yy) * (W{N}.'.*u.values);
            
            % Update the difforder:
            u.diffOrder = u.diffOrder - 1;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  VOLT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = volt(K, u, oneVar)
            %VOLT  Volterra operator with kernel K.
            %   VOLT(K, U) computes the action of the Volterra operator with
            %   kernel K on the Chebyshev interpolant to the points in the
            %   vector U.
            
            persistent x C
            
            % Extract the data:
            N = length(u.values);
            c = diff(u.domain)/2; % Interval scaling.
            
            % Compute cumsum matrix:
            if ( size(C, 2) ~= N )
                C = chebcolloc2.cumsummat(N);
                x = chebpts(N, u.domain);
            end
            
            % Evaluate the kernel:
            if ( nargin == 3 && oneVar )
                KER = K(x);
            else
                [X, Y] = ndgrid(x);
                KER = K(X, Y);
            end
            
            % The Fredholm operator:
            u.values = c * (KER .* C) * (u.values);
            
            % Update the difforder:
            u.diffOrder = u.diffOrder - 1;
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FEVAL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function out = feval(u, y)
            %FEVAL   Evaluate polynomial interpolant of data {X_cheb, U} at a
            % point y using barycentric interpolation.
            
            persistent dom x v
            
            N = length(u.values);
                        
            if ( length(x) ~= N || isempty(dom) || ~all(dom == u.domain) )
                [x, w, v] = chebpts(N, u.domain);
                dom = u.domain;
            end
            
            out = bary(y, u.values, x, v);
            out = reshape(out, size(y));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MISC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Standard Matlab methods. Most of these proceed by simply calling the
        % corresponding method on the values property of the CHEBDOUBLE object.
        
        function u = abs(u)
            u.values = abs(u.values);
        end
        
        function u = acos(u)
            u.values = acos(u.values);
        end
        
        function u = acosd(u)
            u.values = acosd(u.values);
        end
        
        function u = acosh(u)
            u.values = acosh(u.values);
        end
        
        function u = acot(u)
            u.values = acot(u.values);
        end
        
        function u = acotd(u)
            u.values = acotd(u.values);
        end
        
        function u = acoth(u)
            u.values = acoth(u.values);
        end
        
        function u = acsc(u)
            u.values = acsc(u.values);
        end
        
        function u = acscd(u)
            u.values = acscd(u.values);
        end
        
        function u = acsch(u)
            u.values = acsch(u.values);
        end
        
        function u = airy(k, z, varargin)
            if ( nargin == 1 )
                u = k;
                u.values = airy(u.values);
            else
                u = z;
                u.values = airy(k, u.values, varargin{:});
            end
        end
        
        function u = asec(u)
            u.values = asec(u.values);
        end
        
        function u = asecd(u)
            u.values = asecd(u.values);
        end
        
        function u = asech(u)
            u.values = asech(u.values);
        end
        function u = asin(u)
            u.values = asin(u.values);
        end
        
        function u = asind(u)
            u.values = asind(u.values);
        end
        
        function u = asinh(u)
            u.values = asinh(u.values);
        end
        
        function u = atan(u)
            u.values = atan(u.values);
        end
        
        function u = atand(u)
            u.values = atand(u.values);
        end 
        
        function u = atanh(u)
            u.values = atanh(u.values);
        end
        
        function u = besselh(nu, k, u, varargin)
            if ( nargin == 2 )
                u = k;
                u.values = besselh(nu, u.values);
            else
                u.values = besselh(nu, k, u.values, varargin{:});
            end
        end
        
        function u = besseli(nu, u, varargin)
            u.values = besseli(nu, u.values, varargin{:});
        end
        
        function u = besselj(nu, u, varargin)
            u.values = besselj(nu, u.values, varargin{:});
        end
        
        function u = besselk(nu, u, varargin)
            u.values = besselk(nu, u.values, varargin{:});
        end  
        
        function u = bessely(nu, u, varargin)
            u.values = bessely(nu, u.values, varargin{:});
        end 
        
        function u = conj(u)
            u.values = conj(u.values);
        end
        
        function u = cos(u)
            u.values = cos(u.values);
        end
        
        function u = cosd(u)
            u.values = cosd(u.values);
        end
        function u = cosh(u)
            u.values = cosh(u.values);
        end
        
        function u = cot(u)
            u.values = cot(u.values);
        end
        
        function u = cotd(u)
            u.values = cotd(u.values);
        end
        
        function u = coth(u)
            u.values = coth(u.values);
        end
        
        function u = csc(u)
            u.values = csc(u.values);
        end
        
        function u = cscd(u)
            u.values = cscd(u.values);
        end
        
        function u = csch(u)
            u.values = csch(u.values);
        end
        
        function u = ctranspose(u)
            u.values = ctranspose(u.values);
        end
        
        function d = double(u)
            d = zeros(length(u), numel(u));
            for k = 1:numel(u)
                d(:,k) = u(k).values;
            end
        end
        
        function varargout = ellipj(u, M, varargin)
            [varargout{1:nargout}] = ellipj(u.values, M, varargin{:});
            for k = 1:nargout
                u.values = varargout{k};
                varargout{k} = u;
            end
        end
        
        function varargout = ellipke(M, varargin)
            [varargout{1:nargout}] = ellipke(M.values, varargin{:});
            for k = 1:nargout
                M.values = varargout{k};
                varargout{k} = M;
            end
        end 
        
        function u = erf(u)
            u.values = erf(u.values);
        end
        
        function u = erfc(u)
            u.values = erfc(u.values);
        end
        
        function u = erfcinv(u)
            u.values = erfcinv(u.values);
        end
        
        function u = erfcx(u)
            u.values = erfcx(u.values);
        end
        
        function u = erfinv(u)
            u.values = erfinv(u.values);
        end
        
        function u = exp(u)
            u.values = exp(u.values);
        end
        
        function u = expm1(u)
            u.values = expm1(u.values);
        end
        
        function u = extractColumns(u, k)
            u.values = u.values(:, k);
        end
        
        function val = get(u, prop)
            switch prop
                case 'diffOrder'
                    for k = 1:numel(u)
                        val(1,k) = u(k).diffOrder;
                    end
                case 'domain'
                    val = u(1).domain; 
                case 'values'
                    val = u.values;
            end
        end  
        
        function u = heaviside(u)
            u.values = heaviside(u.values);
        end
        
        function u = imag(u)
            u.values = imag(u.values);
        end
        
        function out = isnan(u)
            out = isnan(u.values);
        end
        
        function out = length(u)
            out = length(u(1).values);
        end
        
        function out = mean(f)
            out = sum(f)/diff(f.domain);
        end
        
        function u = minus(u, v)
            u = plus(u, -v);
        end
        
        function u = mrdivide(u, v)
             if ( isnumeric(v) )
                u.values = u.values/v;
            elseif ( isnumeric(u) )
                v.values = u/v.values;
                u = v;
            else
                error('CHEBFUN:CHEBDOUBLE:mrdivide:dim', ...
                    'Matrix dimensions must agree.');
            end
        end
        
        function u = mtimes(u, v)
            if ( isnumeric(v) )
                u.values = u.values*v;
            elseif ( isnumeric(u) )
                v.values = u*v.values;
                u = v;
            elseif ( numel(u.values) < 2 || numel(v.values) < 2 )
                u.values = u.values*v;
                u.diffOrder = max(u.diffOrder, v.diffOrder);   
            else
                error('CHEBFUN:CHEBDOUBLE:mtimes:dim', ...
                    'Matrix dimensions must agree.');
            end
        end
        
        function out = norm(u)
            out = sqrt(sum(u.*u));
        end
        
        function u = rdivide(u, v)
            if ( isnumeric(v) )
                u.values = u.values./v;
            elseif ( isnumeric(u) )
                v.values = u./v.values;
                u = v;
            else
                u.values = u.values./v.values;
                u.diffOrder = max(u.diffOrder, v.diffOrder);
            end
        end
        
        function u = plus(u, v)
            if ( isnumeric(v) )
                u.values = u.values + v;
            elseif ( isnumeric(u) )
                v.values = u + v.values;
                u = v;
            else
                u.values = u.values + v.values;
                u.diffOrder = max(u.diffOrder, v.diffOrder);
            end
        end
        
        function u = power(u, b)
            if ( isnumeric(b) )
                u.values = u.values.^b;
            elseif ( isnumeric(u) )
                b.values = u.*b.values;
                u = b;
            else
                u.values = u.values.^b.values;
                u.diffOrder = max(u.diffOrder, b.diffOrder);
            end
        end
        
        function u = real(u)
            u.values = real(u.values);
        end
        
        function u = sin(u)
            u.values = sin(u.values);
        end
        
        function u = sinc(u)
            u.values = sin(pi*u.values)./(pi*u.values);
            u.values(isnan(u.values)) = 1;
        end
        
        function u = sind(u)
            u.values = sind(u.values);
        end
        
        function u = sinh(u)
            u.values = sinh(u.values);
        end
        
        function u = sqrt(u)
            u.values = sqrt(u.values);
        end
        
        function u = subsref(u, s)
            if ( isnumeric(s.subs{1}) )
                u = feval(u, s.subs{1});
            else
                u = subsref(u.values, s);
            end
        end
        
        function u = tan(u)
            u.values = tan(u.values);
        end
        
        function u = tand(u)
            u.values = tand(u.values);
        end
        
        function u = tanh(u)
            u.values = tanh(u.values);
        end
        
        function u = times(u, v)
            if ( isnumeric(v) )
                u.values = u.values.*v;
            elseif ( isnumeric(u) )
                v.values = u.*v.values;
                u = v;
            elseif ( ~isa(u, 'chebdouble') )
                u = times(feval(u, chebpts(length(v), v.domain)), v);
            elseif ( ~isa(v, 'chebdouble') )
                u = times(u, feval(v, chebpts(length(u), u.domain)));
            else
                u.values = u.values.*v.values;
                u.diffOrder = max(u.diffOrder, v.diffOrder);                
            end
        end
        
        function u = transpose(u)
            u.values = transpose(u.values);
        end
        
        function u = uminus(u)
            u.values = -u.values;
        end
        
        function u = uplus(u)
        end
                
    end
    
end


