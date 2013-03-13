classdef funcheb %< smoothfun % (Abstract)
%FUNCHEB   Approximate smooth functions on [-1,1] with Chebyshev interpolants. 
%
%   Class for approximating smooth functions on the interval [-1,1] using
%   function values Chebyshev points and coefficients of the corresponding
%   1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   FUNCHEB.CONSTRUCTOR(OP) constructs a FUNCHEB object from the function handle
%   OP. OP should be vectorised (i.e., accept a vector input) and ouput a vector
%   of the same length. FUNCHEB objects allow for vectorised construction (i.e.,
%   of multi-valued functions), in which case OP should accept a column vector
%   of length N and return a matrix of size NxM.
%
%   FUNCHEB.CONSTRUCTOR(OP, VSCALE, HSCALE) constructs a FUNCHEB with
%   'happiness' relative to the maximum of the given vertical scale VSCALE
%   (which is updated by the infinity norm of the sampled function values of OP
%   during construction), and the fixed horizontal scale HSCALE. If not given,
%   the VSCALE defaults to 0 initially, and HSCALE defaults to 1.
%
%   FUNCHEB.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF) overrides the default behavior
%   with that given by the preference structure PREF. The constructor will also
%   accept inputs of the form FUNCHEB(OP, PREF), but this usage is not advised.
%   Similarly, one can pass FUNCHEB(OP, VSCALE, HSCALE, EPS), which is
%   equivalent to the call FUNCHEB(OP, VSCALE, HSCALE, FUNCHEB.PREF('eps',EPS)).
%
%   The FUNCHEB class supports construction via interpolation at first- and
%   second-kind Chebyshev points with the classes FUNCHEB1 and FUNCHEB2
%   respectively. The default procedure is to use 2nd-kind points, but this can
%   be overwritted with the preferences PREF = FUNCHEB('tech', cheb1').
%
%   FUNCHEB.CONSTRUCTOR(VALUES, VSCALE, HSCALE, PREF) returns a FUNCHEB object
%   which interpolates the data in the columns of VALUES on a Chebyshev grid.
%   Whether this grid is of first- or second-kind points is determined by
%   PREF.FUNCHEB.tech, as above. FUNCHEB.CONSTRUCTOR({VALUES, COEFFS} allows for
%   the corresponding Chebyshev coefficients to be passed also, and if VALUES is
%   empty the FUNCHEB is constructed directly from the COEFFS.
%
% Examples:
%   % Basic construction:
%   f = funcheb.constructor(@(x) sin(x))
%
%   % Construction with preferences:
%   p = funcheb.pref('tech', 'cheb2'); % See help('funcheb.pref') for details.
%   f = funcheb.constructor(@(x) cos(x), p)
%
%   % Vector-valued construction:
%   f = funcheb.constructor(@(x) [sin(x), cos(x), exp(x)])
%
% See also FUNCHEB.PREF, HAPPINESSCHECK, FUNCHEB1, FUNCHEB2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCHEB Class Description:
%
% The FUNCHEB class is an abstract class for representations of smooth functions
% on the interval [-1,1] via interpolated function values at Chebyshev points
% and coefficients of the corresponding first-kind Chebyshev series expansion.
%
% There are two main instances on the FUNCHEB class; FUNCHEB1 and FUNCHEB2,
% which interpolate on Chebyshev grids of the 1st and 2nd kind respectively.
% Note that although they use different Chebyshev grids in 'value' space, their
% coefficients are always from an expansion in first-kind Chebyshev polynomials
% (i.e., those usualy denoted by $T_k(x)$).
%
% The decision to use FUNCHEB1 or FUNCHEB2 is decided by the funcheb.pref.tech
% property, which should be either of the strings 'cheb1' or 'cheb2'.
%
% The vertical scale VSCALE is used to enforce scale invariance in FUNCHEB
% construction and subsequent operations. For example, that 
%  funcheb.constructor(@(x) 2^300*f(x)) = 2^300*funcheb2.constructor(@(x) f(x)). 
% VSCALE may be optionally passed during to the constructor (if not, it defaults
% to 0), and during construction it is updated to be the maximum magnitude of
% the sampled function values. Similarly the horizontal scale HSCALE is used to
% enforce scale invariance when the input OP has been implicitly mapped from a
% domain other than [-1 1] before being passed to a FUNCHEB constructor.
%
% EPSLEVEL is the happiness level to which the FUNCHEB was constructed (See
% HAPPINESSCHECK.m for full documentation) or a rough accuracy estimate of
% subsequent operations, both relative to VSCALE. Therefore EPSLEVEL could be
% regarded as the number of correct digits in the sampled value that created
% VSCALE. 
%
% Here is a rough guide to how scale and accuracy information is propogated in
% subsequent operations after construction:
%   h = f + c:      
%     h.vscale = max(h.values, [], 1);
%     h.epslevel = (f.epslevel*f.vscale + c*eps)/h.vscale;   % eps(c)/c?
%                  
%   h = f * c:      
%     h.vscale = max(abs(h.values), [], 1) = abs(c)*f.vscale;
%     h.epslevel = f.epslevel + eps; % eps(c)/c?
% 
%   h = f + g:      
%     h.vscale = max(abs(h.values), [], 1);
%     h.epslevel = (f.epslevel*f.vscale + g.epslevel*g.vscale)/h.vscale
% 
%   h = f .* g:
%     h.vscale = max(abs(h.values), [], 1);
%     h.epslevel = (f.epslevel + g.epslevel) * (f.vscale*g.vscale)/h.vscale
% 
%   h = diff(f):
%     h.vscale = max(abs(h.values), [], 1);
%     h.epslevel = n*log(n)f.epslevel*f.vscale; % *(h.vscale/h.vscale)
%     % We don't divide by h.vscale here as we must also multiply by it.
%
% If the input operator OP evaluates to NaN or Inf at any of the sample points
% used by the constructor, then a suitable replacement is found by extrapolating
% (globally) from the numeric values (see EXTRAPOLATE.M). If the preference
% funcheb.pref('extrapolate', TRUE) is set, then the endpoint values -1 and +1
% are always extrapolated (i.e., regardless of whether they evaluate to NaN).
%
% The FUNCHEB classes support the representation of vector-valued functions (for
% example, f = funcheb.constructor(@(x) [sin(x), cos(x)])). In such cases, the
% values and coefficients are stored in a matrix (column-wise), and as such each
% component of the multi-valued function is truncated to the same length, even
% if the demands of 'happiness' imply that one of the components could be
% truncated to a shorter length than the others. All FUNCHEB methods should
% accept such vectorised forms. Note that this representation is distinct from
% an array of FUNCHEB objects, for which there is little to no support.
%
% Class diagram: [<<smoothfun>>] <-- [<<FUNCHEB>>] <-- [funcheb1]
%                                                  <-- [funcheb2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    %% Properties of FUNCHEB objects.
    properties ( Access = public )
        
        % Values (stored in order from left to right) of FUNCHEB at Chebyshev
        % points. The particular Chebyshev points used depend on the instance of
        % the concrete class (1st kind for FUNCHEB1 and 2nd kind for FUNCHEB2).
        % For vectorised FUNCHEB objects, each column represents the
        % interpolated values of a single function.
        values % (nxm double)
        
        % Coefficients in 1st-kind Chebyshev series expansion of the FUNCHEB on
        % [-1,1]. The coefficients are stored in descending order so that c_N is
        % the first entry and c_0 is the last. For vectorised FUNCHEB objects,
        % each column represents the coefficients of a single function.
        coeffs % (nxm double)
        
        % Vertical scale of the FUNCHEB. The is the magnitude of the largest
        % entry in VALUES. It is convenient to store this as a property.
        vscale = 0 % (1xm double >= 0)
        
        % Horizontal scale of the FUNCHEB. Although FUNCHEB objects have in
        % principle no notion of horizontal scale invariance (since they always
        % live on [-1,1]), the input OP has been implicitly mapped. HSCALE is
        % then used to enforce horizontal scale invariance in construction and
        % other subsequent operations that require it. It defaults to 1, and is
        % never updated.
        hscale = 1 % (scalar > 0)        
        
        % Boolean value designating whether the FUNCHEB is 'happy' or not. See
        % HAPPINESSCHECK.m for full documentation.
        ishappy % (logical)
        
        % Happiness level to which the FUNCHEB was constructed (See
        % HAPPINESSCHECK.m for full documentation) or a rough accuracy estimate
        % of subsequent operations (See FUNCHEB class documentaion for details).
        epslevel % (double >= 0)
    end
    
    %% CLASS CONSTRUCTOR:
    methods (Static)
        function obj = constructor(op, vscale, hscale, pref)
            % Constructor for the FUNCHEB class.
            
            % We can't return an empty FUNCHEB, so pass an empty OP down.
            if ( nargin == 0  )
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
            
            % Obtain preferences.
            if ( nargin == 2 && isstruct(vscale) )
                % vscale was actually a preference.
                pref = funcheb.pref(vscale);
                vscale = 0;
                hscale = 1;
            elseif ( nargin == 3 && isstruct(hscale) )
                % hscale was actually a preference.
                pref = funcheb.pref(hscale);
                hscale = 1;
            elseif ( nargin < 4 )
                % Create:
                pref = funcheb.pref;
            elseif ( ~isstruct(pref) )
                % An eps was passed.
                pref = funcheb.pref('eps', pref);
            else
                % Merge:
                pref = funcheb.pref(pref);
            end

            % Call the relevent constructor:
            if ( strcmpi(pref.funcheb.tech, 'cheb1') )
                % Construct:
                obj = funcheb1(op, vscale, hscale, pref);
            else
                % Construct:
                obj = funcheb2(op, vscale, hscale, pref);
            end
            
        end
    end
    
    
    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY THIS CLASS. 
    methods (Abstract)
       
        % Compose method. (Not implemented here as refinement is defined also).
        h = compose(f, op, g, pref)

        % Get method. [TODO]: Requirement should be inherited from smoothfun.
        val = get(f, prop);
        
        % Set method. [TODO]: Requirement should be inherited from smoothfun.
%         f = set(f, prop, val); % [TODO]: Do we actually need a set method?
        
    end
    
    %% ABSTRACT STATIC METHODS REQUIRED BY THIS CLASS.
    methods (Abstract, Static)
        
        % Compute Chebyshev barycentric weights.
        w = barywts(n)
        
        % Convert values to coefficients.
        coeffs = chebpoly(values)
        
        % Convert coefficients to values.
        values = chebpolyval(coeffs)
        
        % Compute Chebyshev points (x) and optionally quadrature (w) and
        % barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % Extrapolate (for NaNs / endpoints).
        [values, maskNaN, maskInf] = extrapolate(values)
        
        % Make a FUNCHEB. (Constructor shortcut)
        f = make(varargin);
        
        % Refinement function for FUNCHEB construction. (Evaluates OP on grid)
        [values, opints, giveUp] = refine(op, values, pref)
        
        % Compute Chebyshev quadrature weights.
        w = quadwts(n)
        
        % Test this class.
        pass = test(); % [TODO]: This should be inherited from smoothfun.
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods 
        
        % Convert an array of FUNCHEB objects into a vector-valued FUNCHEB.
        f = cell2mat(f)
        
        % Plot (semilogy) the Chebyshev coefficients of a FUNCHEB object.
        h = chebpolyplot(f, varargin)
        
        % Check the happiness of a FUNCHEB. (Classic definition).
        [ishappy, epslevel, cutoff] = classicCheck(f, pref)
        
        % Complex conjugate of a FUNCHEB.
        f = conj(f)

        % Indefinite integral of a FUNCHEB.
        f = cumsum(f, pref)

        % Derivative of a FUNCHEB.
        f = diff(f, k, dim)

        % Evaluate a FUNCHEB.
        y = feval(f,x)
        
        % Flip columns of a vectorised FUNCHEB object.
        f = fliplr(f)
        
        % Happiness test for a FUNCHEB
        [ishappy, epslevel, cutoff] = happinessCheck(f, op, pref)
        
        % Imaginary part of a FUNCHEB.
        f = imag(f)

        % Compute the inner product of two FUNCHEB objects.
        out = innerProduct(f, g)

        % True for an empty FUNCHEB.
        out = isempty(f)

        % Test if FUNCHEB objects are equal.
        out = isequal(f, g)

        % Test if a FUNCHEB is bounded.
        out = isfinite(f)
        
        % Test if a FUNCHEB is unbounded.
        out = isinf(f)

        % Test if a FUNCHEB is has any NaN values.
        out = isnan(f)
        
        % True for real FUNCHEB.
        out = isreal(f)

        % Length of a FUNCHEB.
        len = length(f)

        % [TODO]: Implement looseCheck.
        [ishappy, epslevel, cutoff] = looseCheck(f, pref)
        
        % Convert a vector-valued FUNCHEB into an ARRAY of FUNCHEB objects.
        g = mat2cell(f, M, N)
        
        % Global maximum of a FUNCHEB on [-1,1].
        [maxVal, maxPos] = max(f)
        
        % Global minimum of a FUNCHEB on [-1,1].
        [minVal, minPos] = min(f)
        
        % Global minimum and maximum on [-1,1].
        [vals, pos] = minandmax(f)

        % Subtraction of two FUNCHEB objects.
        f = minus(f, g)

        % Left matrix divide for FUNCHEB objects.
        X = mldivide(A, B)

        % Right matrix divide for a FUNCHEB.
        X = mrdivide(B, A)
        
        % Multiplication of FUNCHEB objects.
        f = mtimes(f, c)
        
        % Basic linear plot for FUNCHEB objects.
        varargout = plot(f, varargin)
        
        % Addition of two FUNCHEB objects.
        f = plus(f, g)
        
        % Return the points used by a FUNCHEB.
        out = points(f)

        % Polynomial coefficients of a FUNCHEB.
        out = poly(f)

        % Populate a FUNCHEB class with values.
        f = populate(f, op, vscale, hscale, pref)
        
        % QR factorisation of a multivalued FUNCHEB.
        [f, R, E] = qr(f, flag)
        
        % Right array divide for a FUNCHEB.
        f = rdivide(f, c, pref)

        % Real part of a FUNCHEB.
        f = real(f)

        % Restrict a FUNCHEB to a subinterval.
        f = restrict(f, s)

        % Roots of a FUNCHEB in the interval [-1,1].
        out = roots(f, varargin)
        
        % Test an evaluation of the input OP against a FUNCHEB approx.
        pass = sampleTest(op, f)
        
        % Trim trailing Chebyshev coefficients of a FUNCHEB object. 
        f = simplify(f, pref, force)
        
        % Size of a FUNCHEB.
        [siz1, siz2] = size(f, varargin)

        % [TODO]: Implement strictCheck
        [ishappy, epslevel, cutoff] = strictCheck(f, pref)
        
        % Definite integral of a FUNCHEB on the interval [-1,1].
        out = sum(f, dim)
        
        % FUNCHEB multiplication.
        f = times(f, g, varargin)
        
        % Unary minus of a FUNCHEB.
        f = uminus(f)

        % Unary plus of a FUNCHEB.
        f = uplus(f)

    end
    
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS.
    methods ( Static = true ) 
          
        % Evaluate a Chebyshev polynomial using barycentric interpolation.
        fx = bary(x, gvals, xk, vk, kind, p)

        % Clenshaw's algorithm for evaluating a Chebyshev polynomial.
        out = clenshaw(x, coeffs)  
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        
    end
    
end


