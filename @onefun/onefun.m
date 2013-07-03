classdef onefun % (Abstract) 
%ONEFUN     Approximate smooth functions on [-1,1]. 
%   Abstract (interface) class for approximating functions on the interval
%   [-1,1].
%
% Constructor inputs:
%   ONEFUN.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF) creates a representation of the
%   operator OP defined on the interval [-1,1]. If PREF.ONEFUN.BLOWUP = 1 (which
%   is the default) then the ONEFUN constructor calls
%    SINGFUN.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF2), else it calls
%    SMOOTHFUN.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF2) if PREF.ONEFUN.BLOWUP = 0, 
%   where in both cases PREF2 is a merge of PREF with the default SINGFUN and
%   ONEFUN preferences, respectively.
%
%   See SMOOTHFUN and SINGFUN for further documentation.
%
% See also ONEFUN.pref, SINGFUN, SMOOTHFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONEFUN Class Description:
%
% The ONEFUN class is an abstract class for representations of functions on the
% interval [-1,1].
%
% The current instances of ONEFUNs are SMOOTHFUNs and SINGFUNs. The former are
% used to represenet smooth functions on [-1,1}, whereas the latter are able to
% represent some forms of endpoint singularites. 
%
% Class diagram: [<<fun>>] <>-- [<<ONEFUN>>] <-- [<<smoothfun>>]
%                                            <-- [singfun]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    %% Properties of ONEFUN objects.
    properties ( Access = public )

        % Vertical scale of the ONEFUN. This is an estimate for the magnitude of
        % the largest entry of the ONEFUN. For an array-valued ONEFUN, VSCALE is
        % an row vector and the kth entry corresponds to the kth column in the
        % ONEFUN.
        vscale = 0 % (1xm double >= 0)

        % Horizontal scale of the ONEFUN. Although ONEFUN objects have in
        % principle no notion of horizontal scale invariance (since they always
        % live on [-1,1]), the input OP may have been implicitly mapped. HSCALE
        % is then used to enforce horizontal scale invariance in construction
        % and other subsequent operations that require it. It defaults to 1 and
        % is never updated.
        hscale = 1 % (scalar > 0)

        % Boolean value designating whether the ONEFUN is a 'happy'
        % representation or not. (See subclasses for further documentation.)
        ishappy % (logical)

        % Happiness level to which the ONEFUN was constructed, or a rough
        % accuracy estimate of subsequent operations. (See subclasses for
        % further documentation.)
        epslevel % (double >= 0)
    end

    methods (Static)
        function obj = constructor(op, vscale, hscale, pref)
            
            % We can't return an empty ONEFUN, so pass an empty OP down.
            if ( nargin == 0 )
                op = [];
            end
            
            % Define vscale if none given:
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            % Define hscale if none given:
            if ( nargin < 3 || isempty(hscale) )
                hscale = 1;
            end
            % Determine preferences if not given, merge if some are given:
            if ( nargin < 4 || isempty(pref) )
                pref = onefun.pref;
            else
                pref = onefun.pref(pref);
            end
       
            % Call the relevent constructor:
            if ( isa(op, 'onefun') )
                % OP is already a ONEFUN!
                obj = op;
                
            elseif ( pref.onefun.blowup )
                % BLOWUP mode; Call SINGFUN.
                
                % Merge preferences:
                pref = singfun.pref(pref, pref.onefun);
                
                % Call singfun constructor:
                obj = singfun(op, vscale, hscale, pref);
                
                % Return just a FUN if no singularities found:
                if ( ~any(obj.exps) )
                    obj = obj.smoothfun; 
                end 
                
            else
                % STANDARD mode; Call SMOOTHFUN.
                
                % Merge preferences:
                pref = smoothfun.pref(pref, pref.onefun);
                
                % Call SMOOTHFUN constructor:
                obj = smoothfun.constructor(op, vscale, hscale, pref);
                
            end
        
        end
        
    end
    
    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY ONEFUN CLASS.
    methods ( Abstract = true )
        
        % Convert an array of ONEFUN objects into a array-valued ONEFUN.
        f = cell2mat(f)

        % Complex conjugate of a ONEFUN.
        f = conj(f)
        
        % ONEFUN obects are not transposable.
        f = ctranspose(f)

        % Indefinite integral of a ONEFUN.
        f = cumsum(f, m, pref)

        % Derivative of a ONEFUN.
        f = diff(f, k, dim)
        
        % Evaluate a ONEFUN.
        y = feval(f, x)

        % Flip columns of an array-valued ONEFUN object.
        f = fliplr(f)
        
        % Flip/reverse a ONEFUN object.
        f = flipud(f)

        % Imaginary part of a ONEFUN.
        f = imag(f)

        % Compute the inner product of two ONEFUN objects.
        out = innerProduct(f, g)

        % True for an empty ONEFUN.
        out = isempty(f)

        % Test if ONEFUN objects are equal.
        out = isequal(f, g)

        % Test if a ONEFUN is bounded.
        out = isfinite(f)

        % Test if a ONEFUN is unbounded.
        out = isinf(f)

        % Test if a ONEFUN has any NaN values.
        out = isnan(f)

        % True for real ONEFUN.
        out = isreal(f)
        
        % True for zero ONEFUN objects
        out = iszero(f)
        
        % Length of a ONEFUN.
        len = length(f)

        % Convert a array-valued ONEFUN into an ARRAY of ONEFUN objects.
        g = mat2cell(f, M, N)

        % Global maximum of a ONEFUN on [-1,1].
        [maxVal, maxPos] = max(f)

        % Global minimum of a ONEFUN on [-1,1].
        [minVal, minPos] = min(f)

        % Global minimum and maximum on [-1,1].
        [vals, pos] = minandmax(f)

        % Subtraction of two ONEFUN objects.
        f = minus(f, g)

        % Left matrix divide for ONEFUN objects.
        X = mldivide(A, B)

        % Right matrix divide for a ONEFUN.
        X = mrdivide(B, A)

        % Multiplication of ONEFUN objects.
        f = mtimes(f, c)
        
        % Compute a Legendre series expansion of a ONEFUN object:
        c = legpoly(f)

        % Basic linear plot for ONEFUN objects.
        varargout = plot(f, varargin)
        
        % Obtain data used for plotting a ONEFUN object:
        data = plotData(f)

        % Addition of two ONEFUN objects.
        f = plus(f, g)

        % Polynomial coefficients of a ONEFUN.
        out = poly(f)

        % QR factorisation of an array-valued ONEFUN.
        [f, R, E] = qr(f, flag, methodFlag)

        % Right array divide for a ONEFUN.
        f = rdivide(f, c, pref)

        % Real part of a ONEFUN.
        f = real(f)

        % Restrict a ONEFUN to a subinterval.
        f = restrict(f, s)

        % Roots of a ONEFUN in the interval [-1,1].
        out = roots(f, varargin)

        % Simplify the representation of a ONEFUN.
        f = simplify(f, pref, force)

        % Size of a ONEFUN.
        [siz1, siz2] = size(f, varargin)

        % Definite integral of a ONEFUN on the interval [-1,1].
        out = sum(f, dim)

        % ONEFUN multiplication.
        f = times(f, g, varargin)
        
        % ONEFUN obects are not transposable.
        f = transpose(f)

        % Unary minus of a ONEFUN.
        f = uminus(f)

        % Unary plus of a ONEFUN.
        f = uplus(f)

    end

    %% ABSTRACT STATIC METHODS REQUIRED BY ONEFUN CLASS.
    methods ( Abstract = true, Static = true )
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
    end
    
    %% Methods implimented by ONEFUN class.
    methods 
        
    end
    
    %% Static methods implimented by ONEFUN class.
    methods ( Static = true ) 
        
    end
    
end
