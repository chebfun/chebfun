classdef bndfun < fun
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
%   that given by the preference structure PREF. See CHEBPREF for details.
%
%   BNDFUN(VALUES, DOMAIN, VSCALE, HSCALE, PREF) returns a BNDFUN object with a
%   ONEFUN constructed by the data in the columns of VALUES (if supported by
%   ONEFUN class constructor).
%
% See ONEFUN for further documentation of the ONEFUN class.
%
% See also FUN, CHEBPREF, ONEFUN.

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
% Class diagram: [<<FUN>>] <>-- [<<onefun>>]
%                   ^
%                   |  
%                [bndfun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% CLASS CONSTRUCTOR:
    methods
        
        function obj = bndfun(op, domain, vscale, hscale, pref)
            
            % Construct an empty bndfun:
            if ( (nargin == 0) || isempty(op) )
                return
            end
            
            % Obtain preferences if none given:
            if ( (nargin < 5) || isempty(pref))
                pref = chebpref();
            else
                pref = chebpref(pref);
            end
            
            % Use default domain if none given. Otherwise, check whether the
            % domain input has correct dimensions
            if ( (nargin < 2) || isempty(domain) )
                domain = pref.domain;
            elseif ( ~all(size(domain) == [1, 2]) ) || diff(domain) <= 0
                error('CHEBFUN:BNDFUN:domain',...
                    ['Domain argument should be a row vector with two ', ...
                    'entries in increasing order.']);
            end
            
            % Check domain:
            if ( any(isinf(domain)) )
                error('CHEBFUN:BNDFUN:UNBND',...
                    'Should not encounter unbounded domain in bndfun class.');
            elseif ( ~all(size(domain) == [1 2]) )
                error('CHEBFUN:BNDFUN:UNBND',...
                    'Domain should be a 1x2 vector.');
            end
            
            % Define scales if none given:
            if ( (nargin < 3) || isempty(vscale) )
                vscale = 0;
            end

            if ( (nargin < 4) || isempty(hscale) )
                % [TODO]: Or should this be 1? What does the chebfun level pass
                % down?
                hscale = norm(domain, inf); 
            end

            linmap = bndfun.createMap(domain);
            % Include linear mapping from [-1,1] to [a,b] in the op:
            if ( isa(op, 'function_handle') && ~all(domain == [-1, 1]) && ...
                    ~isnumeric(op) )
                op = @(x) op(linmap.for(x));
            end
            
            % Call the ONEFUN constructor:
            obj.onefun = onefun.constructor(op, vscale, hscale/diff(domain), pref);
            
            % Add the domain and mapping:
            obj.domain = domain;
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

        % Compose a BNDFUN with an operator or another BNDFUN.
        f = compose(f, op, g, pref)
        
        % Indefinite integral of a BNDFUN.
        [f, rVal] = cumsum(f, m, dim, shift)
        
        % Derivative of a BNDFUN.
        f = diff(f, k, dim)
       
        % Change of domains of BNDFUN via linear change of variables.
        f = changeMap(f,newdom)
        
        % Evaluate a BNDFUN.
        y = feval(f, x)
        
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