classdef chebtech < smoothfun % (Abstract)
%CHEBTECH   Approximate smooth functions on [-1,1] with Chebyshev interpolants.
%
%   Abstract class for approximating smooth functions on the interval [-1,1]
%   using function values at Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% See also CHEBTECH1, CHEBTECH2, CHEBTECH.TECHPREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBTECH Class Description:
%
% The CHEBTECH class is an abstract class for representations of smooth
% functions on the interval [-1,1] via interpolated function values at Chebyshev
% points and coefficients of the corresponding first-kind Chebyshev series
% expansion.
%
% There are two concrete realizations of the CHEBTECH class--CHEBTECH1 and
% CHEBTECH2--which interpolate on Chebyshev grids of the 1st and 2nd kind,
% respectively.  Note that although they use different Chebyshev grids in
% 'value' space, their coefficients are always from an expansion in first-kind
% Chebyshev polynomials (i.e., those usually denoted by $T_k(x)$).
%
% The vertical scale VSCALE is used to enforce scale invariance in CHEBTECH
% construction and subsequent operations. For example, that
%   chebtech1(@(x) 2^300*f(x)) = 2^300*chebtech1(@(x) f(x)).
%
% VSCALE may be optionally passed to the constructor (if not, it defaults to 0),
% and during construction it is updated to be the maximum magnitude of the
% sampled function values. Similarly the horizontal scale HSCALE is used to
% enforce scale invariance when the input OP has been implicitly mapped from a
% domain other than [-1 1] before being passed to a CHEBTECH constructor.
%
% If the input operator OP in a call to a concrete CHEBTECH constructor, say,
% CHEBTECH1(OP), evaluates to NaN or Inf at any of the sample points used by the
% constructor, then a suitable replacement is found by extrapolating (globally)
% from the numeric values (see EXTRAPOLATE.M). If the EXTRAPOLATE preference is
% set to TRUE (See CHEBTECH.TECHPREF), then the endpoint values -1 and +1 are
% always extrapolated (i.e., regardless of whether they evaluate to NaN).
%
% The CHEBTECH classes support the representation of array-valued functions (for
% example, f = chebtech1(@(x) [sin(x), cos(x)])). In such cases, the values and
% coefficients are stored in a matrix (column-wise), and as such each component
% of the array-valued function is truncated to the same length, even if the
% demands of 'happiness' imply that one of the components could be truncated to
% a shorter length than the others. All CHEBTECH methods should accept such
% array-valued forms. Note that this representation is distinct from an array of
% CHEBTECH objects, for which there is little to no support.
%
% Class diagram: [<<SMOOTHFUN>>] <-- [<<CHEBTECH>>] <-- [CHEBTECH1]
%                                                   <-- [CHEBTECH2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )

        % Coefficients in 1st-kind Chebyshev series expansion of the CHEBTECH on
        % [-1,1]. The coefficients are stored in descending order so that c_N is
        % the first entry and c_0 is the last. For array-valued CHEBTECH
        % objects, each column represents the coefficients of a single function.
        coeffs % (nxm double)

        % Boolean value designating whether the CHEBTECH is 'happy' or not. See
        % HAPPINESSCHECK.m for full documentation.
        ishappy % (logical)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false, Abstract = true )

        % Compose method. (Not implemented here as refinement is defined also).
        h = compose(f, op, g, data, pref)

        % Get method.
        val = get(f, prop);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true, Abstract = true )
        
        % Alias Chebyshev coefficients.
        coeffs = alias(coeffs, m)
        
        % Angles of Chebyshev points. (i.e., acos(chebpts(n))
        t = angles(n)

        % Compute Chebyshev barycentric weights.
        w = barywts(n)

        % Compute Chebyshev points (x) and optionally quadrature (w) and
        % barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % Convert coefficients to values.
        values = coeffs2vals(coeffs)

        % Make a CHEBTECH. (Constructor shortcut)
        f = make(varargin);

        % Refinement function for CHEBTECH construction. (Evaluates OP on grid)
        [values, points, giveUp] = refine(op, values, pref)

        % Compute Chebyshev quadrature weights.
        w = quadwts(n)
        
        % Convert values to coefficients.
        coeffs = vals2coeffs(values)

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS (IMPLEMENTED BY THIS CLASS.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )

        % Clenshaw's algorithm for evaluating a Chebyshev polynomial.
        out = clenshaw(x, coeffs)

        % Convert Chebyshev-T coefficients to Chebyshev-U coefficients.
        cU = chebTcoeffs2chebUcoeffs(cT)

        % Retrieve and modify preferences for this class.
        p = techPref(q)

        % Parse inputs passed to the constructor via the DATA argument.
        data = parseDataInputs(data, pref)
    end

end
