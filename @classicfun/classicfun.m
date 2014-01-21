classdef classicfun < fun % (Abstract)
%CLASSICFUN   Represent global functions on an interval [a, b].
%
%   Abstract (interface) class for representing global functions on an interval
%   [a, b], which can either be bounded or unbounded. Functions are
%   approximated via a ONEFUN object, which lives on the interval [-1, 1],
%   stored in the CLASSICFUN.  Forward and inverse maps stored in the CLASSICFUN object map
%   the interval [-1, 1] to [a, b], and vice versa.
%
% Constructor inputs:
%   CLASSICFUN.CONSTRUCTOR(OP, DOMAIN) constructs a CLASSICFUN object from the function handle
%   OP by mapping the DOMAIN to [-1, 1], and constructing an ONEFUN object to
%   represent the function prescribed by OP. DOMAIN should be a row vector with
%   two elements in increasing order. OP should be vectorised (i.e., accept a
%   vector input) and output a vector of the same length as its input.
%
%   CLASSICFUN.CONSTRUCTOR(OP, DOMAIN, VSCALE, HSCALE) allows the constructor of the
%   ONEFUN of the CLASSICFUN to use information about vertical and horizontal scales.
%   If not given (or given as empty), the VSCALE defaults to 0 initially, and
%   HSCALE defaults to 1.
%
%   CLASSICFUN.CONSTRUCTOR(OP, DOMAIN, VSCALE, HSCALE, PREF) overrides the default
%   behavior with that given by the CHEBPREF object PREF. See CHEBPREF for
%   details.
%
%   CLASSICFUN.CONSTRUCTOR(VALUES, DOMAIN, VSCALE, HSCALE, PREF) returns a CLASSICFUN object
%   with a ONEFUN constructed by the data in the columns of VALUES (if supported
%   by ONEFUN class constructor).
%
% See ONEFUN for further documentation of the ONEFUN class.
%
% See also CHEBPREF, ONEFUN, BNDFUN, UNBNDFUN.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLASSICFUN Class Description:
%
% The CLASSICFUN class is an abstract class for representations of functions on the
% interval [a, b]. It achieves this by taking a ONEFUN on [-1, 1] and applying
% a mapping.
%
% The current instances of CLASSICFUNs are BNDFUNs and UNBNDFUNs. The former are used
% to represent functions on bounded domains, whereas the latter are able to
% represent some functions on unbounded domains.
%
% Note that all binary CLASSICFUN operators (methods which can take two CLASSICFUN arguments)
% assume that the domains of the CLASSICFUN objects agree. The methods will not throw
% warnings in case the domains don't agree, but their output will not be
% meaningful.
%
% Class diagram: [chebfun] <>-- [<<FUN>>] <>-- [<<onefun>>]
%                                 ^   ^
%                                /     \
%                          [bndfun]   [unbndfun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Properties of CLASSICFUN objects.
    properties (Access = public)
        % The domain of the CLASSICFUN object, [a, b].
        domain

        % The mapping which maps [-1, 1] to [a, b], and vice versa.
        mapping

        % The ONEFUN of the CLASSICFUN object, which does the actual approximation of
        % the function the CLASSICFUN object represents. The ONEFUN object lives on the
        % interval [-1, 1], the mapping of the CLASSICFUN object takes care of mapping
        % between [-1, 1] and vice versa.
        onefun
    end
    
    %% CLASS CONSTRUCTOR:
    methods (Static = true)
        function obj = constructor(op, domain, vscale, hscale, pref)
            
            % We can't return an empty CLASSICFUN, so pass an empty OP down.
            if ( nargin == 0  )
                op = [];
            end
            
            % Obtain preferences if none given:
            if ( nargin < 5 )
                pref = chebpref();
            else
                pref = chebpref(pref);
            end
            
            % Get domain if none given:
            if ( nargin < 2 || isempty(domain) )
                domain = pref.domain;
            end
            
            % Get vscale if none given:
            if ( nargin < 3 || isstruct(vscale) )
                vscale = 0;
            end
            
            % Get hscale if none given:
            if ( nargin < 4 || isempty(vscale) )
                hscale = norm(domain, inf);
            end

            % [TODO]: Explain this. Only becomes relevant with UNBNDFUN
            if ( isinf(hscale) )
                hscale = 1;
            end

            % Call constructor depending on domain:
            if ( ~any(isinf(domain)) )
                % Construct a BNDFUN object:
                obj = bndfun(op, domain, vscale, hscale, pref);
                
            else
                % Construct an UNBNDFUN object:
                obj = unbndfun(op, domain, vscale, hscale, pref);
                
            end
            
        end

    end
    
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS.
    methods (Static = true)

    end
    
    %% ABSTRACT STATIC METHODS REQUIRED BY THIS CLASS.
    methods(Abstract = true, Static = true)
                
        % Map from [-1, 1] to the domain of the CLASSICFUN.
        m = createMap(domain);  
        
        % Make a CLASSICFUN. (Constructor shortcut)
        f = make(varargin);
    end
    
    %% ABSTRACT METHODS REQUIRED BY THIS CLASS.
    methods(Abstract = true)
        % [TODO]: Once UNBNDFUN and CHEBFUN advance, we should revisit this
        % list, and add/throw away abstract methods as appropriate.
        
        % Compose a CLASSICFUN with an operator or another CLASSICFUN
        f = compose(f, op, g, pref)
        
        % Indefinite integral of a CLASSICFUN.
        f = cumsum(f, m, pref)
        
        % Derivative of a CLASSICFUN.
        f = diff(f, k, dim)
        
        % Evaluate a CLASSICFUN.
        y = feval(f, x)
        
        % Compute the inner product of two CLASSICFUN objects.
        out = innerProduct(f, g)
        
        % Data for plotting a CLASSICFUN
        data = plotData(f, g);
        
        % Restrict a CLASSICFUN to a subinterval.
        f = restrict(f, s)
        
        % Definite integral of a CLASSICFUN on its domain of definition.
        out = sum(f, dim)
        
    end           
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % Extract information for DISPLAY.
        info = dispData(f)
        
        % Extract columns of an array-valued FUN object.
        f = extractColumns(f, columnIndex);

        % Round a FUN towards zero.
        g = fix(f);
        
        % Round a FUN towards minus infinity.
        g = floor(f);

        % Flip columns of an array-valued FUN object.
        f = fliplr(f)
        
        % Get properties of a FUN.
        out = get(f, prop);
        
        % Imaginary part of a FUN.
        f = imag(f)

        % True for an empty FUN.
        out = isempty(f)

        % Test if FUN objects are equal.
        out = isequal(f, g)

        % Test if a FUN is bounded.
        out = isfinite(f)

        % Test if a FUN is unbounded.
        out = isinf(f)

        % Test if a FUN has any NaN values.
        out = isnan(f)

        % True for real FUN.
        out = isreal(f)
        
        % Test if a FUN object is built upon SINGFUN.
        out = issing(f)
        
        % True for zero FUN objects
        out = iszero(f)
        
        % Length of a FUN.
        len = length(f)

        % FUN logical.
        f = logical(f)

        % Convert an array-valued FUN into a cell array of FUN objects.
        g = mat2cell(f, M, N)

        % Global maximum of a FUN on [a,b].
        [maxVal, maxPos] = max(f)

        % Global minimum of a FUN on [a,b].
        [minVal, minPos] = min(f)

        % Global minimum and maximum on [a,b].
        [vals, pos] = minandmax(f)

        % Subtraction of two FUN objects.
        f = minus(f, g)

        % Multiplication of FUN objects.
        f = mtimes(f, c)

        % FUN logical NOT.
        f = not(f)

        % FUN logical OR.
        h = or(f, g)

        % Basic linear plot for FUN objects.
        varargout = plot(f, varargin)
        
        % 3-D plot for FUN objects.
        varargout = plot3(f, g, h, varargin)

        % Addition of two FUN objects.
        f = plus(f, g)

        % Right array divide for a FUN.
        f = rdivide(f, c, pref)

        % Real part of a FUN.
        f = real(f)

        % Roots of a FUN in the interval [a,b].
        out = roots(f, varargin)
        
        % Round a FUN towards nearest integer.
        g = round(f)
        
        % Signum of a FUN. (f should have no zeros in its domain)
        f = sign(f, pref)

        % Simplify the ONEFUN of a FUN object.
        f = simplify(f, tol)

        % Size of a FUN.
        [size1, size2] = size(f, varargin)

        % FUN multiplication.
        f = times(f, g, varargin)
        
        % FUN objects are not transposable.
        f = transpose(f)

        % Unary minus of a FUN.
        f = uminus(f)

        % Unary plus of a FUN.
        f = uplus(f)

    end
end
