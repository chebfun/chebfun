classdef classicfun < fun % (Abstract)
%CLASSICFUN   Represent global functions on an interval [a, b].
%
%   Abstract (interface) class for representing global functions on an interval
%   [a, b], which can either be bounded or unbounded. Functions are
%   approximated via a ONEFUN object, which lives on the interval [-1, 1],
%   stored in the CLASSICFUN.  Forward and inverse maps stored in the
%   CLASSICFUN object map the interval [-1, 1] to [a, b], and vice versa.
%
% Constructor inputs:
%   CLASSICFUN.CONSTRUCTOR(OP) constructs a CLASSICFUN object from the function
%   handle OP on the domain determined by the default stored in CHEBFUNPREF by
%   mapping the domain to [-1, 1], and constructing an ONEFUN object to
%   represent the function prescribed by OP.
%
%   CLASSICFUN.CONSTRUCTOR(OP, DATA) does the same but uses the domain
%   specified by DATA.DOMAIN.  If DATA or DATA.DOMAIN is empty or if DATA has
%   no DOMAIN field, the default domain from CHEBFUNPREF is used.  The DATA
%   structure is passed as-is to the appropriate concrete subclass constructor.
%
%   CLASSICFUN.CONSTRUCTOR(OP, DATA, PREF) overrides the default behavior with
%   that given by the preferences in the structure or CHEBFUNPREF object PREF.
%   See CHEBFUNPREF for details.
%
% See ONEFUN for further documentation of the ONEFUN class.
%
% See also CHEBFUNPREF, ONEFUN, BNDFUN, UNBNDFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLASSICFUN Class Description:
%
% The CLASSICFUN class is an abstract class for representations of functions on
% the interval [a, b]. It achieves this by taking a ONEFUN on [-1, 1] and
% applying a mapping.
%
% The current instances of CLASSICFUNs are BNDFUNs and UNBNDFUNs. The former
% are used to represent functions on bounded domains, whereas the latter are
% able to represent some functions on unbounded domains.
%
% Note that all binary CLASSICFUN operators (methods which can take two
% CLASSICFUN arguments) assume that the domains of the CLASSICFUN objects
% agree. The methods will not throw warnings in case the domains don't agree,
% but their output will not be meaningful.
%
% Class diagram: [CHEBFUN]<>---[<<FUN>>]------[<<CLASSICFUN>>] <>-- [<<onefun>>]
%                                  |            /   ^     ^
%                             [DELTAFUN]<>-----/   /       \
%                                               [bndfun]   [unbndfun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % The domain of the CLASSICFUN object, [a, b].
        domain

        % The mapping which maps [-1, 1] to [a, b], and vice versa.
        mapping

        % The ONEFUN of the CLASSICFUN object, which does the actual
        % approximation of the function the CLASSICFUN object represents. The
        % ONEFUN object lives on the interval [-1, 1], the mapping of the
        % CLASSICFUN object takes care of mapping between [-1, 1] and vice
        % versa.
        onefun
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        function obj = constructor(op, data, pref)
            
            % Parse inputs.
            if ( nargin < 1  )
                % We can't return an empty CLASSICFUN, so pass an empty OP down.
                op = [];
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

            % Call constructor depending on domain:
            if ( ~any(isinf(data.domain)) )
                % Construct a BNDFUN object:
                obj = bndfun(op, data, pref);
            else
                % Construct an UNBNDFUN object:
                obj = unbndfun(op, data, pref);
            end
            
        end
        
    end
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT METHODS REQUIRED BY THIS CLASS.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Abstract = true )
        % [TODO]: Once UNBNDFUN and CHEBFUN advance, we should revisit this
        % list, and add/throw away abstract methods as appropriate.
        
        % Compose a CLASSICFUN with an operator or another CLASSICFUN
        f = compose(f, op, g, data, pref)
        
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC ABSTRACT METHODS REQUIRED BY THIS CLASS.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true, Abstract = true )
                
        % Map from [-1, 1] to the domain of the CLASSICFUN.
        m = createMap(domain);  
        
        % Make a CLASSICFUN. (Constructor shortcut)
        f = make(varargin);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE METHODS (IMPLEMENTED BY THIS ABSTRACT CLASS.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )

        function out = isPeriodicTech(f)
        %ISPERIODICTECH   Test if the smooth part of f is is constructed with a
        %basis of periodic functions.
            out = isPeriodicTech(f.onefun);
        end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS IMPLEMENTED IN THIS FILE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = parseDataInputs(data, pref)
%PARSEDATAINPUTS   Parse inputs from the DATA structure and assign defaults.

if ( ~isfield(data, 'domain') || isempty(data.domain) )
    data.domain = pref.domain;
end

end
