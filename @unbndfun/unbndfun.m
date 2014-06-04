classdef unbndfun < classicfun
%UNBNDFUN    Represent global functions on an unbounded interval [-inf inf] or
% a semi-infinite interval [-inf b] or [a inf].
%
% Constructor inputs:
%   UNBNDFUN(OP, DOMAIN) constructs an UNBNDFUN object from the function handle
%   OP by mapping the DOMAIN to [-1 1], and constructing an ONEFUN object to
%   represent the function prescribed by OP. DOMAIN should be a row vector with
%   two elements in increasing order and at least one entry of this two-entry
%   vector should be inf or -inf. OP should be vectorised (i.e., accept a vector
%   input) and output a vector of the same length as its input.
%
%   UNBNDFUN(OP, DOMAIN, VSCALE, HSCALE) allows the constructor of the ONEFUN of
%   the UNBNDFUN to use information about vertical and horizontal scales. If not
%   given (or given as empty), the VSCALE defaults to 0 initially, and HSCALE
%   defaults to 1.
%
%   UNBNDFUN(OP, DOMAIN, VSCALE, HSCALE, PREF) overrides the default behavior
%   with that given by the preference structure PREF. See CHEBFUNPREF.
%
%   UNBNDFUN(VALUES, DOMAIN, VSCALE, HSCALE, PREF) returns a UNBNDFUN object
%   with a ONEFUN constructed by the data in the columns of VALUES (if supported
%   by ONEFUN class constructor).
%
% See ONEFUN for further documentation of the ONEFUN class.
%
% See also FUN, CHEBFUNPREF, ONEFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNBNDFUN Class Description:
%
% The UNBNDFUN class represents global functions on the infinite or
% semi-infinite interval [-inf, b], [a, inf], or [-inf, inf]. It achieves this
% by taking a onefun on [-1, 1] and applying a nonlinear mapping.
%
% Note that all binary UNBNDFUN operators (methods which can take two UNBNDFUN
% arguments) assume that the domains of the UNBNDFUN objects agree. The methods
% will not throw warnings if assumption is violated, but the results will not be
% meaningful under that circumstance.
%
% Class diagram: [<<CLASSICFUN>>] <>-- [<<onefun>>]
%                    ^
%                    |
%                [unbndfun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% CLASS CONSTRUCTOR:   
    methods
        
        function obj = unbndfun(op, data, pref)
            % Parse inputs.
            if ( (nargin < 1) || isempty(op) )
                obj.domain = [];
                obj.mapping = [];
                obj.onefun = [];
                return
            end

            if ( (nargin < 2) || isempty(data) )
                    data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = chebfunpref();
            else
                pref = chebfunpref(pref);
            end

            data = parseDataInputs(data, pref);

            % Check domain
            if ( ~all(size(data.domain) == [1, 2]) || (diff(data.domain) <= 0) )
                error('CHEBFUN:UNBNDFUN:domain',...
                    ['Domain argument should be a row vector with two ', ...
                    'entries in increasing order.']);
            elseif ( ~any(isinf(data.domain)) )
                error('CHEBFUN:UNBNDFUN:boundedDomain',...
                    'Domain argument should be unbounded.');
            end

            % Remap the OP to be a function on [-1, 1].
            unbndmap = unbndfun.createMap(data.domain);
            if ( isa(op, 'function_handle') )
                op = @(x) op(unbndmap.for(x));
            elseif ( isnumeric(op) )
                if ( ~any(op(:)) )
                    op = @(x) zeros(length(x), size(op, 2));
                else
                    %[TODO]: Implement this.
                    error('CHEBFUN:UNBNDFUN:inputValues',...
                        ['UNBNDFUN does not support non-zero construction ' ...
                         'from values.']);
                end
            end

            % Deal with exponents for singular functions.
            if ( isempty(data.exponents) )
                % The user hasn't supplied exponents, but functions on
                % unbounded domains are often singular once remapped to [-1,
                % 1].  Try to detect this by seeing if the function is infinite
                % at either endpoint and, if so, enable singularity detection.
                lVal = feval(op, -1);
                rVal = feval(op, 1);
                if ( any(isinf([lVal rVal])) )
                    pref.enableSingularityDetection = true;
                    singType = pref.singPrefs.defaultSingType;
                    data.singType = {singType, singType};
                end
            else
                % Remapping to [-1, 1] negates exponents, which are given with
                % respect to the function on the infinite domain.
                ind = isinf(data.domain);
                pref.enableSingularityDetection = true;
                data.exponents(ind) = -data.exponents(ind);
            end

            % Call the ONEFUN constructor:
            obj.onefun = onefun.constructor(op, data, pref);

            % Add the domain and mapping:
            obj.domain = data.domain;
            obj.mapping = unbndmap;
        end
    end
    
    %% STATIC METHODS IMPLEMENTED BY UNBNDFUN CLASS.
    methods ( Static = true ) 
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin);
        
        % Noninear map from [-1, 1] to the domain of the UNBNDFUN.
        m = createMap(domain);
        
        % Make a UNBNDFUN (constructor shortcut):
        f = make(varargin); 
        
    end
       
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % Compose an UNBNDFUN with an operator or another FUN.
        f = compose(f, op, g, pref)
        
        % Indefinite integral of an UNBNDFUN.
        [f, rVal] = cumsum(f, dim)
        
        % Derivative of an UNBNDFUN.
        f = diff(f, k, dim)
       
        % Change of domains of an UNBNDFUN via linear change of variables.
        f = changeMap(f,newdom)
        
        % Evaluate an UNBNDFUN.
        y = feval(f, x, varargin)
        
        % Flip/reverse an UNBNDFUN object.
        f = flipud(f)
        
        % Compute the inner product of two UNBNDFUN objects.
        out = innerProduct(f, g)
        
        % Left matrix divide for UNBNDFUN objects.
        X = mldivide(A, B)

        % Right matrix divide for an UNBNDFUN.
        X = mrdivide(B, A)
        
        % Estimate the Inf-norm of an UNBNDFUN
        out = normest(f);
        
        % Data for plotting an UNBNDFUN
        data = plotData(f, g);
                
        % Polynomial coefficients of an UNBNDFUN.
        out = poly(f)
        
        % UNBNDFUN power function.
        f = power(f, b)
        
        % QR factorisation for UNBNDFUN object is not supported.
        [f, R, E] = qr(f, flag)

        % Restrict an UNBNDFUN to a subinterval.
        f = restrict(f, s)
        
        % Roots of an UNBNDFUN in an unbounded domain.
        out = roots(f, varargin)
        
        % Definite integral of an UNBNDFUN on the its domain.
        out = sum(f, dim)
    end    
end

function data = parseDataInputs(data, pref)

if ( ~isfield(data, 'domain') || isempty(data.domain) )
    data.domain = pref.domain;
end

if ( ~isfield(data, 'hscale') || isempty(data.hscale) )
    % TODO:  Why is the hscale of an unbounded domain always 1?
    data.hscale = 1;
end

if ( ~isfield(data, 'exponents') || isempty(data.exponents) )
    data.exponents = [];
end

if ( ~isfield(data, 'singType') || isempty(data.singType) )
    data.singType = [];
end

end
