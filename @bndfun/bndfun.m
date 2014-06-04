classdef bndfun < classicfun
%BNDFUN   Represent global functions on a bounded interval [a, b].
%
%   Class for representing global functions on a bounded interval [a, b].
%   Functions are approximated via a ONEFUN object, which lives on the interval
%   [-1, 1], stored in the BNDFUN. Forward and inverse maps stored in the
%   BNDFUN object map the interval [-1, 1] to [a, b], and vice versa.
%
% Constructor inputs:
%   BNDFUN(OP, DOMAIN) constructs a BNDFUN object from the function handle OP by
%   mapping the DOMAIN to [-1, 1], and constructing an ONEFUN object to
%   represent the function prescribed by OP. DOMAIN should be a row vector with
%   two elements in increasing order. OP should be vectorised (i.e., accept a
%   vector input) and output a vector of the same length as its input.
%
%   BNDFUN(OP, DOMAIN, VSCALE, HSCALE) allows the constructor of the ONEFUN of
%   the BNDFUN to use information about vertical and horizontal scales. If not
%   given (or given as empty), the VSCALE defaults to 0 initially, and HSCALE
%   defaults to 1.
%
%   BNDFUN(OP, DOMAIN, VSCALE, HSCALE, PREF) overrides the default behavior with
%   that given by the preference structure PREF. See CHEBFUNPREF for details.
%
%   BNDFUN(VALUES, DOMAIN, VSCALE, HSCALE, PREF) returns a BNDFUN object with a
%   ONEFUN constructed by the data in the columns of VALUES (if supported by
%   ONEFUN class constructor).
%
% See ONEFUN for further documentation of the ONEFUN class.
%
% See also CLASSICFUN, CHEBFUNPREF, ONEFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BNDFUN Class Description:
%
% The BNDFUN class is a class for representations of globally smooth functions
% on a finite interval [a, b]. It achieves this by taking a ONEFUN on [-1, 1]
% and applying a linear mapping to re-scale its domain.
%
% Note that all binary BNDFUN operators (methods which can take two BNDFUN
% arguments) assume that the domains of the BNDFUN objects agree. The methods
% will not issue warnings if this condition is violated, but the results will
% not be meaningful.
%
% Class diagram: [<<CLASSICFUN>>] <>-- [<<onefun>>]
%                   ^
%                   |  
%                [bndfun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% CLASS CONSTRUCTOR:
    methods
        function obj = bndfun(op, data, pref)
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

            % Check the domain input.
            if ( ~all(size(data.domain) == [1, 2]) || (diff(data.domain) <= 0) )
                error('CHEBFUN:BNDFUN:badDomain', ...
                    ['Domain argument should be a row vector with two ', ...
                    'entries in increasing order.']);
            elseif ( any(isinf(data.domain)) )
                error('CHEBFUN:BNDFUN:unboundedDomain', ...
                    'Should not encounter unbounded domain in bndfun class.');
            end

            % TODO:  Why do we rescale the hscale like this?
            data.hscale = data.hscale / diff(data.domain);

            % Remap the OP to be a function on [-1, 1].
            linmap = bndfun.createMap(data.domain);
            if ( isa(op, 'function_handle') && ~all(data.domain == [-1, 1]) ...
                    && ~isnumeric(op) )
                op = @(x) op(linmap.for(x));
            end

            % Call the ONEFUN constructor:
            obj.onefun = onefun.constructor(op, data, pref);

            % Add the domain and mapping:
            obj.domain = data.domain;
            obj.mapping = linmap;
        end
    end
    
    %% STATIC METHODS IMPLEMENTED BY BNDFUN CLASS.
    methods ( Static = true ) 

        % Linear map from [-1, 1] to the domain of the BNDFUN.
        m = createMap(domain);
        
        % Make a BNDFUN (constructor shortcut):
        f = make(varargin);
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods
        
        % Convolution of BNDFUN F with BNDFUN G.
        h = conv(f, g)

        % Compose a BNDFUN with an operator or another BNDFUN.
        f = compose(f, op, g, pref)
        
        % Indefinite integral of a BNDFUN.
        [f, rVal] = cumsum(f, dim)
        
        % Derivative of a BNDFUN.
        f = diff(f, k, dim)
       
        % Change of domains of BNDFUN via linear change of variables.
        f = changeMap(f,newdom)
        
        % Evaluate a BNDFUN.
        y = feval(f, x, varargin)
        
        % Flip/reverse a BNDFUN object.
        f = flipud(f)

        % Compute the inner product of two BNDFUN objects.
        out = innerProduct(f, g)
        
        % Left matrix divide for BNDFUN objects.
        X = mldivide(A, B)

        % Right matrix divide for a BNDFUN.
        X = mrdivide(B, A)
        
        % Data for plotting a BNDFUN
        data = plotData(f, g, h);
                
        % Polynomial coefficients of a BNDFUN.
        out = poly(f)
        
        % BNDFUN power function.
        f = power(f, b)
        
        % QR factorisation of an array-valued BNDFUN.
        [f, R, E] = qr(f, flag)
        
        % Restrict a BNDFUN to a subinterval.
        f = restrict(f, s)
        
        % Definite integral of a BNDFUN on the interval [a, b].
        out = sum(f, dim)
    end
end

function data = parseDataInputs(data, pref)

if ( ~isfield(data, 'domain') || isempty(data.domain) )
    data.domain = pref.domain;
end

if ( ~isfield(data, 'hscale') || isempty(data.hscale) )
    % TODO:  Or should this be 1?  What does chebfun pass down?
    data.hscale = norm(data.domain, inf);
end

end
