classdef chebtech < smoothfun % (Abstract)
%CHEBTECH   Approximate smooth functions on [-1,1] with Chebyshev interpolants.
%
%   Class for approximating smooth functions on the interval [-1,1] using
%   function values at Chebyshev points and coefficients of the corresponding
%   1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   CHEBTECH.CONSTRUCTOR(OP) constructs a CHEBTECH object from the function
%   handle OP by evaluating it on an increasingly fine set of Chebyshev points
%   (see below). OP should be vectorised (i.e., accept a vector input) and
%   output a vector of the same length. CHEBTECH objects allow for array-valued
%   functions, in which case OP should accept a column vector of length N and
%   return a matrix of size NxM.
%
%   CHEBTECH.CONSTRUCTOR(OP, VSCALE, HSCALE) constructs a CHEBTECH with
%   'happiness' relative to the maximum of the given vertical scale VSCALE
%   (which is updated by the infinity norm of the sampled function values of OP
%   during construction), and the fixed horizontal scale HSCALE. If not given
%   (or given as empty), the VSCALE defaults to 0 initially, and HSCALE defaults
%   to 1.
%
%   CHEBTECH.CONSTRUCTOR(OP, VSCALE, HSCALE, PREF) overrides the default behavior
%   with that given by the preference structure PREF. See CHEBTECH.pref for
%   details. The CHEBTECH class supports construction via interpolation at first-
%   and second-kind Chebyshev points with the classes CHEBTECH1 and CHEBTECH2
%   respectively. The default procedure is to use 2nd-kind points, but this can
%   be overwritten with the preferences PREF = CHEBTECH('tech', 'cheb1').
%
%   CHEBTECH.CONSTRUCTOR(VALUES, VSCALE, HSCALE, PREF) returns a CHEBTECH object
%   which interpolates the data in the columns of VALUES on a Chebyshev grid.
%   Whether this grid is of first- or second-kind points is determined by
%   PREF.CHEBTECH.TECH, as above. CHEBTECH.CONSTRUCTOR({VALUES, COEFFS}, ...)
%   allows for the corresponding Chebyshev coefficients to be passed also, and
%   if VALUES is empty the CHEBTECH is constructed directly from the COEFFS. No
%   adaptivity takes place with this form of construction, but VALUES are still
%   checked for happiness. If COEFFS are passed, the resulting CHEBTECH is
%   always deemed 'happy'.
%
% Examples:
%   % Basic construction:
%   f = chebtech.constructor(@(x) sin(x))
%
%   % Construction with preferences:
%   p = chebtech.pref('tech', 'cheb2'); % See HELP('chebtech.pref') for details.
%   f = chebtech.constructor(@(x) cos(x), [], [], p)
%
%   % Array-valued construction:
%   f = chebtech.constructor(@(x) [sin(x), cos(x), exp(x)])
%
% See also CHEBTECH.PREF, HAPPINESSCHECK, CHEBTECH1, CHEBTECH2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBTECH Class Description:
%
% The CHEBTECH class is an abstract class for representations of smooth functions
% on the interval [-1,1] via interpolated function values at Chebyshev points
% and coefficients of the corresponding first-kind Chebyshev series expansion.
%
% There are two concrete realizations of the CHEBTECH class--CHEBTECH1 and
% CHEBTECH2--which interpolate on Chebyshev grids of the 1st and 2nd kind,
% respectively.  Note that although they use different Chebyshev grids in
% 'value' space, their coefficients are always from an expansion in first-kind
% Chebyshev polynomials (i.e., those usually denoted by $T_k(x)$).
%
% The decision to use CHEBTECH1 or CHEBTECH2 is decided by the CHEBTECH.PREF.TECH
% property, which should be either of the strings 'cheb1' or 'cheb2'.
%
% The vertical scale VSCALE is used to enforce scale invariance in CHEBTECH
% construction and subsequent operations. For example, that
%
% chebtech.constructor(@(x) 2^300*f(x)) = 2^300*chebtech2.constructor(@(x) f(x)).
%
% VSCALE may be optionally passed to the constructor (if not, it defaults to
% 0), and during construction it is updated to be the maximum magnitude of the
% sampled function values. Similarly the horizontal scale HSCALE is used to
% enforce scale invariance when the input OP has been implicitly mapped from a
% domain other than [-1 1] before being passed to a CHEBTECH constructor.
%
% EPSLEVEL is the happiness level to which the CHEBTECH was constructed (See
% HAPPINESSCHECK.m for full documentation) or a rough accuracy estimate of
% subsequent operations, both relative to VSCALE. Therefore EPSLEVEL could be
% regarded as the number of correct digits in the sampled value that created
% VSCALE.
%
% Here is a rough guide to how scale and accuracy information is propagated in
% subsequent operations after construction:
%   h = f + c:
%     h.vscale = max(h.values, [], 1);
%     h.epslevel = (f.epslevel*f.vscale + eps(c))/h.vscale;
%
%   h = f * c:
%     h.vscale = max(abs(h.values), [], 1) = abs(c)*f.vscale;
%     h.epslevel = f.epslevel + eps(c)/c;
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
%     % [TODO]: Figure this out rigourously.
%     h.epslevel = n*log(n)f.epslevel*f.vscale; % *(h.vscale/h.vscale)
%     % We don't divide by h.vscale here as we must also multiply by it.
%
%   h = cumsum(f):
%     h.vscale = max(abs(h.values), [], 1);
%     [TODO]: h.epslevel = ???
%
% If the input operator OP evaluates to NaN or Inf at any of the sample points
% used by the constructor, then a suitable replacement is found by extrapolating
% (globally) from the numeric values (see EXTRAPOLATE.M). If the preference
% CHEBTECH.PREF('extrapolate', TRUE) is set, then the endpoint values -1 and +1
% are always extrapolated (i.e., regardless of whether they evaluate to NaN).
%
% The CHEBTECH classes support the representation of array-valued functions (for
% example, f = chebtech.constructor(@(x) [sin(x), cos(x)])). In such cases, the
% values and coefficients are stored in a matrix (column-wise), and as such each
% component of the array-valued function is truncated to the same length, even
% if the demands of 'happiness' imply that one of the components could be
% truncated to a shorter length than the others. All CHEBTECH methods should
% accept such array-valued forms. Note that this representation is distinct from
% an array of CHEBTECH objects, for which there is little to no support.
%
% Class diagram: [<<smoothfun>>] <-- [<<CHEBTECH>>] <-- [chebtech1]
%                                                  <-- [chebtech2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Properties of CHEBTECH objects.
    properties ( Access = public )

        % Values of CHEBTECH at Chebyshev points (stored in order from left to
        % right). The particular Chebyshev points used depend on the instance
        % of the concrete class (1st kind for CHEBTECH1 and 2nd kind for
        % CHEBTECH2).  % For array-valued CHEBTECH objects, each column
        % represents the interpolated values of a single function.
        values % (nxm double)

        % Coefficients in 1st-kind Chebyshev series expansion of the CHEBTECH on
        % [-1,1]. The coefficients are stored in descending order so that c_N is
        % the first entry and c_0 is the last. For array-valued CHEBTECH
        % objects, each column represents the coefficients of a single function.
        coeffs % (nxm double)

        % Vertical scale of the CHEBTECH. This is a row vector storing the
        % magnitude of the largest entry in each column of VALUES. It is
        % convenient to store this as a property.
%         vscale = 0 % (1xm double >= 0)

        % Horizontal scale of the CHEBTECH. Although CHEBTECH objects have in
        % principle no notion of horizontal scale invariance (since they always
        % live on [-1,1]), the input OP may have been implicitly mapped.
        % HSCALE is then used to enforce horizontal scale invariance in
        % construction and other subsequent operations that require it. It
        % defaults to 1 and is never updated.
%         hscale = 1 % (scalar > 0)

        % Boolean value designating whether the CHEBTECH is 'happy' or not. See
        % HAPPINESSCHECK.m for full documentation.
%         ishappy % (logical)

        % Happiness level to which the CHEBTECH was constructed (See
        % HAPPINESSCHECK.m for full documentation) or a rough accuracy estimate
        % of subsequent operations (See CHEBTECH class documentation for details).
%         epslevel % (double >= 0)
    end

    %% CLASS CONSTRUCTOR:
    methods ( Static = true )
        function obj = constructor(op, vscale, hscale, pref)
            % Constructor for the CHEBTECH class.

            % We can't return an empty CHEBTECH, so pass an empty OP down.
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
            % Determine preferences if not given, merge if some are given:
            if ( nargin < 4 || isempty(pref) )
                pref = chebtech.pref;
            else
                pref = chebtech.pref(pref);
            end

            % Call the relevant constructor:
            if ( strcmpi(pref.chebtech.tech, 'cheb1') )
                % Construct:
                obj = chebtech1(op, vscale, hscale, pref);
            else
                % Construct:
                obj = chebtech2(op, vscale, hscale, pref);
            end

        end
    end


    %% ABSTRACT (NON-STATIC) METHODS REQUIRED BY THIS CLASS.
    methods ( Abstract = true )

        % Compose method. (Not implemented here as refinement is defined also).
        h = compose(f, op, g, pref)

        % Get method.
        val = get(f, prop);

    end

    %% ABSTRACT STATIC METHODS REQUIRED BY THIS CLASS.
    methods ( Abstract = true, Static = true )
        
        % Alias Chebyshev coefficients.
        coeffs = alias(coeffs, m)

        % Compute Chebyshev barycentric weights.
        w = barywts(n)

        % Convert values to coefficients.
        coeffs = chebpoly(values)

        % Convert coefficients to values.
        values = chebpolyval(coeffs)

        % Compute Chebyshev points (x) and optionally quadrature (w) and
        % barycentric (v) weights.
        [x, w, v] = chebpts(n)

        % Make a CHEBTECH. (Constructor shortcut)
        f = make(varargin);

        % Refinement function for CHEBTECH construction. (Evaluates OP on grid)
        [values, points, giveUp] = refine(op, values, pref)

        % Compute Chebyshev quadrature weights.
        w = quadwts(n)

    end

    %% METHODS IMPLEMENTED BY THIS CLASS.
    methods

        % Convert an array of CHEBTECH objects into a array-valued CHEBTECH.
        f = cell2mat(f)

        % Plot (semilogy) the Chebyshev coefficients of a CHEBTECH object.
        h = chebpolyplot(f, varargin)

        % Check the happiness of a CHEBTECH. (Classic definition).
        [ishappy, epslevel, cutoff] = classicCheck(f, pref)

        % Complex conjugate of a CHEBTECH.
        f = conj(f)
        
        % CHEBTECH obects are not transposable.
        f = ctranspose(f)

        % Indefinite integral of a CHEBTECH.
        f = cumsum(f, m, pref)

        % Derivative of a CHEBTECH.
        f = diff(f, k, dim)
        
        % Extrapolate (for NaNs / Infs).
        [values, maskNaN, maskInf] = extrapolate(f)

        % Evaluate a CHEBTECH.
        y = feval(f, x)

        % Flip columns of an array-valued CHEBTECH object.
        f = fliplr(f)
        
        % Flip/reverse a CHEBTECH object.
        f = flipud(f)

        % Happiness test for a CHEBTECH
        [ishappy, epslevel, cutoff] = happinessCheck(f, op, pref)

        % Imaginary part of a CHEBTECH.
        f = imag(f)

        % Compute the inner product of two CHEBTECH objects.
        out = innerProduct(f, g)

        % True for an empty CHEBTECH.
        out = isempty(f)

        % Test if CHEBTECH objects are equal.
        out = isequal(f, g)

        % Test if a CHEBTECH is bounded.
        out = isfinite(f)

        % Test if a CHEBTECH is unbounded.
        out = isinf(f)

        % Test if a CHEBTECH has any NaN values.
        out = isnan(f)

        % True for real CHEBTECH.
        out = isreal(f)
        
        % True for zero CHEBTECH objects
        out = iszero(f)
        
        % Length of a CHEBTECH.
        len = length(f)

        % A 'loose' (i.e., not too strict) check for happiness.
        [ishappy, epslevel, cutoff] = looseCheck(f, pref)

        % Convert a array-valued CHEBTECH into an ARRAY of CHEBTECH objects.
        g = mat2cell(f, M, N)

        % Global maximum of a CHEBTECH on [-1,1].
        [maxVal, maxPos] = max(f)

        % Global minimum of a CHEBTECH on [-1,1].
        [minVal, minPos] = min(f)

        % Global minimum and maximum on [-1,1].
        [vals, pos] = minandmax(f)

        % Subtraction of two CHEBTECH objects.
        f = minus(f, g)

        % Left matrix divide for CHEBTECH objects.
        X = mldivide(A, B)

        % Right matrix divide for a CHEBTECH.
        X = mrdivide(B, A)

        % Multiplication of CHEBTECH objects.
        f = mtimes(f, c)
        
        % Compute a Legendre series expansion of a CHEBTECH object:
        c = legpoly(f)

        % Basic linear plot for CHEBTECH objects.
        varargout = plot(f, varargin)
        
        % Obtain data used for plotting a CHEBTECH object:
        data = plotData(f)

        % Addition of two CHEBTECH objects.
        f = plus(f, g)

        % Return the points used by a CHEBTECH.
        out = points(f)

        % Polynomial coefficients of a CHEBTECH.
        out = poly(f)

        % Populate a CHEBTECH class with values.
        f = populate(f, op, vscale, hscale, pref)
        
        % Adjust the number of points used in a CHEBTECH.
        f = prolong(f, n)

        % QR factorisation of an array-valued CHEBTECH.
        [f, R, E] = qr(f, flag, methodFlag)

        % Right array divide for a CHEBTECH.
        f = rdivide(f, c, pref)

        % Real part of a CHEBTECH.
        f = real(f)

        % Restrict a CHEBTECH to a subinterval.
        f = restrict(f, s)

        % Roots of a CHEBTECH in the interval [-1,1].
        out = roots(f, varargin)

        % Test an evaluation of the input OP against a CHEBTECH approx.
        pass = sampleTest(op, f)

        % Trim trailing Chebyshev coefficients of a CHEBTECH object.
        f = simplify(f, pref, force)

        % Size of a CHEBTECH.
        [siz1, siz2] = size(f, varargin)

        % Strict happiness check.
        [ishappy, epslevel, cutoff] = strictCheck(f, pref)

        % Definite integral of a CHEBTECH on the interval [-1,1].
        out = sum(f, dim)

        % CHEBTECH multiplication.
        f = times(f, g, varargin)
        
        % CHEBTECH obects are not transposable.
        f = transpose(f)

        % Unary minus of a CHEBTECH.
        f = uminus(f)

        % Unary plus of a CHEBTECH.
        f = uplus(f)

    end

    %% STATIC METHODS IMPLEMENTED BY THIS CLASS.
    methods ( Static = true )

        % Evaluation using the barycentric interpolation formula.
        fx = bary(x, gvals, xk, vk)

        % Clenshaw's algorithm for evaluating a Chebyshev polynomial.
        out = clenshaw(x, coeffs)

        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)

    end

end
