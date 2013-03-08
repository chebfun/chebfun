classdef funcheb1 < funcheb
%FUNCHEB1 Approximate smooth functions on [-1,1] with Chebyshev interpolants.
%
%   Class for approximating smooth functions on the interval [-1,1]
%   using function values at 1st-kind Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   FUNCHEB1(OP) constructs a FUNCHEB1 object from the function handle OP. OP
%   should be vectorized (i.e., accept a vector input) and ouput a vector of
%   the same length. FUNCHEB1 objects allow for vectorised construction (i.e.,
%   of multi-valued function), in which case OP should accept a vector of length
%   N and return a matrix of size NxM.
%
%   FUNCHEB1(OP, VSCALE) constructs a FUNCHEB1 with 'happiness' (see
%   funcheb1.HAPPINESSCHECK.m) relative to the maximum of the given vertical scale (VSCALE)
%   and the (column-wise) infinity norm of the sampled function values of OP. If not given,
%   the VSCALE defaults to 0 initially.
%
%   FUNCHEB1(OP, VSCALE, HSCALE) uses a 'happiness' to both the vertical scale
%   VSCALE (as above) and the horizontal scale HSCALE. If not given, this
%   defaults to 1.
%
%   FUNCHEB1(OP, VSCALE, HSCALE, PREF) overrides the default behavior with that
%   given by the preference structure PREF. The constructor will also accept
%   inputs of the form FUNCHEB1(OP, PREF), but this usage is not advised.
%   Similarly, one can pass FUNCHEB1(OP, VSCALE, PREF). Furthermore, one can
%   replace PREF by TOL, the desired tolerance of the construction, which is
%   equivelent to passing a PREF with PREF.FUNCHEB1.eps = TOL.
%
%   FUNCHEB1(VALUES, ...) returns a FUNCHEB1 object which interpolates the
%   values in the columns of VALUES at 1st-kind Chebyshev points and
%   FUNCHEB1({VALUES, COEFFS}, ... ) uses the Chebyshev coefficients passed in
%   COEFFS rather than computing them.
%
% Examples:
%   % Basic construction
%   f = funcheb1(@(x) sin(x))
%
%   % Construction with preferences
%   p = funcheb1.pref('sampletest', 0); % See help('funcheb1.pref') for details
%   f = funcheb1(@(x) sin(x), p)
%
%   % Vector-valued construction:
%   f = funcheb1(@(x) [sin(x), cos(x), exp(x)])
%
% See also FUNCHEB1.pref, FUNCHEB1.chebpts, FUNCHEB1.happinesscheck, FUNCHEB1.refine.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCHEB1 Class Description:
%
% The FUNCHEB1 class represents smooth functions on the interval [-1,1] using
% function values at 1st-kind Chebyshev points and coefficients of the
% corresponding first-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 1st-kind grid (see funcheb1.REFINE.m) until the
% representation is deemed 'happy' (see funcheb1.HAPPINESSCHECK.m). The resulting
% object can be used to evaluate and operate on the input function.
%
% The vertical scale VSCALE is used to enforce scale invariance in FUNCHEB1
% construction and subsequent operations. For example, that 
%          funcheb1(@(x) 2^300*f(x)) = 2^300*funcheb1(@(x) f(x)). 
% VSCALE may be optionally passed during to the constructor (if not, it
% defaults to 0), and during construction it is updated to be the maximum
% magnitude of the sampled function values.
%
% If the input operator OP evaluates to NaN or Inf at any of the sample points
% used by the constructor, then a suitable replacement is found by extrapolating
% (globally) from the numeric values (see EXTRAPOLATE.M). If the preference
% funcheb1.pref('extrapolate', TRUE) is set, then the endpoint values -1 and +1
% are always extrapolated (i.e., regardless of whether they evaluate to NaN).
%
% The FUNCHEB1 class supports the representation of vector-valued functions
% (for example, f = funcheb1(@(x) [sin(x), cos(x)])). In such cases, the values
% and coefficients are stored in a matrix (column-wise), and as such each
% component of the multi-valued function is truncated to the same length, even
% if the demands of 'happiness' imply that one of the components could be
% truncated to a shorter length than the others. All FUNCHEB1 methods should
% accept such vectorised forms. Note that this representation is distinct from
% an array of funcheb1 objects, such as [funcheb1(@(x) sin(x), funcheb1(@(x)
% cos(x)], for which there is little to no support.
%
% Class diagram: [funcheb] <-- [funcheb1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Methods implemented by this m-file.
    methods
        
        function obj = funcheb1(op, vscale, hscale, pref)
            %Constructor for the FUNCHEB1 class.
            
            % Return an empty FUNCHEB1 on null input:
            if ( nargin == 0 || isempty(op) )
                return
            end
           
            % Define vscale if none given:
            if ( nargin < 2 || isempty(vscale) )
                vscale = 0;
            end
            if ( nargin < 3 || isempty(hscale) )
                hscale = 1;
            end
            
            % Obtain preferences:
            if ( nargin == 2 && isstruct(vscale) )
                % vscale was actually a preference.
                pref = funcheb1.pref(vscale);
                vscale = 0;
                hscale = 1;
            elseif ( nargin == 3 && isstruct(hscale) )
                % hscale was actually a preference.
                pref = funcheb1.pref(hscale);
                hscale = 1;                
            elseif ( nargin < 4 )
                % Create:
                pref = funcheb1.pref;
            elseif ( ~isstruct(pref) )
                % An eps was passed.
                pref = funcheb1.pref('eps', pref);
            else
                % Merge:
                pref = funcheb1.pref(pref);
            end
            
            % Force nonadaptive construction if pref.funcheb1.n is numeric:
            if ( ~isempty(pref.funcheb1.n) && ~isnan(pref.funcheb1.n) )
                % Evaluate the op on the Chebyshev grid of given size:
                op = feval(op, funcheb1.chebpts(pref.funcheb1.n));
            end
            
            % Actual construction takes place here:
            obj = populate(obj, op, vscale, hscale, pref);
            
            % Check for NaNs: (if not happy)
            if ( ~obj.ishappy )
                % Check for NaNs:
                if ( any(any(isnan(obj.values(:)))) )
                    error('CHEBFUN:FUNCHEB1:constructor:naneval', ...
                        'Function returned NaN when evaluated.')
                end
            elseif ( ~obj.ishappy )
                % Here we throw an error if NaNs were encountered anywhere.
                if ( any(isnan(obj.values(:))) )
                    error('CHEBFUN:FUNCHEB1:constructor:naneval2', ...
                        'Function returned NaN when evaluated.')
                end
            end
            
        end
        
    end
    
    %% Static methods implemented by FUNCHEB1 class.
    % (This list is alphabetical)
    methods ( Static = true )
        
        % aliasing.
        coeffs = alias(coeffs, m)
        
        % Evaluate a Chebyshev interpolant using barycentric formula.
        out = bary(x, values, kind)
        
        % Compute Chebyshev barycentric weights.
        w = barywts(n)

        % Convert values to coefficients.
        coeffs = chebpoly(values)
        
        % Convert coefficients to values.
        values = chebpolyval(coeffs)
        
        % Compute Chebyshev points of first kind (x) and optionally
        % quadrature (w) and barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % Extrapolate (for NaNs / endpoints).
        [values, maskNaN, maskInf] = extrapolate(values)
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Refinement function for FUNCHEB1 construction. (Evaluates OP on grid)
        [values, opints, giveUp] = refine(op, values, pref)
        
        % Compute Chebyshev quadrature weights.
        w = quadwts(n)
        
        % Test the FUNCHEB1 class.
        pass = test(varargin);
        
        % Make a FUNCHEB1. (Constructor shortcut)
        f = make(varargin);

    end
    
    %% Private static methods implemented by FUNCHEB1 class.
    methods ( Access = private, Static = true )
        
    end
    
end
