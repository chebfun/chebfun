classdef chebtech1 < chebtech
%CHEBTECH1   Approximate smooth functions on [-1,1] with Chebyshev interpolants.
%
%   Class for approximating smooth functions on the interval [-1,1]
%   using function values at 1st-kind Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   CHEBTECH1(OP) constructs a CHEBTECH1 object from the function handle OP. OP
%   should be vectorized (i.e., accept a vector input) and ouput a vector of
%   the same length. CHEBTECH1 objects allow for vectorised construction (i.e.,
%   of multi-valued function), in which case OP should accept a vector of length
%   N and return a matrix of size NxM.
%
%   CHEBTECH1(OP, VSCALE) constructs a CHEBTECH1 with 'happiness' (see
%   CHEBTECH1.HAPPINESSCHECK) relative to the maximum of the given vertical scale (VSCALE)
%   and the (column-wise) infinity norm of the sampled function values of OP. If not given,
%   the VSCALE defaults to 0 initially.
%
%   CHEBTECH1(OP, VSCALE, HSCALE) uses a 'happiness' to both the vertical scale
%   VSCALE (as above) and the horizontal scale HSCALE. If not given (or given as
%   empty), this defaults to 1.
%
%   CHEBTECH1(OP, VSCALE, HSCALE, PREF) overrides the default behavior with that
%   given by the preference structure PREF.
%
%   CHEBTECH1(VALUES, ...) returns a CHEBTECH1 object which interpolates the
%   values in the columns of VALUES at 1st-kind Chebyshev points and
%   CHEBTECH1({VALUES, COEFFS}, ... ) uses the Chebyshev coefficients passed in
%   COEFFS rather than computing them.
%
% Examples:
%   % Basic construction:
%   f = chebtech1(@(x) sin(x))
%
%   % Construction with preferences:
%   p = chebtech.pref('sampletest', 0); % See help('chebtech.pref') for details
%   f = chebtech1(@(x) sin(x), p)
%
%   % Vector-valued construction:
%   f = chebtech1(@(x) [sin(x), cos(x), expchebtech.pref(x)])
%
% See also CHEBTECH.PREF, CHEBPTS, HAPPINESSCHECK, REFINE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [TODO]: There is great and probably unnecessary duplication between
%   a lot of this documentation and what you find in CHEBTECH2.  
% CHEBTECH1 Class Description:
%
% The CHEBTECH1 class represents smooth functions on the interval [-1,1] using
% function values at 1st-kind Chebyshev points and coefficients of the
% corresponding 1st-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 1st-kind grid (see CHEBTECH1.REFINE) until the
% representation is deemed 'happy' (see CHEBTECH1.HAPPINESSCHECK). The resulting
% object can be used to evaluate and operate on the input function.
%
% The vertical scale VSCALE is used to enforce scale invariance in CHEBTECH1
% construction and subsequent operations, for example, to enforce that 
%          chebtech1(@(x) 2^300*f(x)) = 2^300*chebtech1(@(x) f(x)). 
%
% VSCALE may be optionally passed to the constructor (if not, it
% defaults to 0), and during construction it is updated to be the maximum
% magnitude of the sampled function values.
%
% If the input operator OP evaluates to NaN or Inf at any of the sample points
% used by the constructor, then a suitable replacement is found by extrapolating
% (globally) from the numeric values (see EXTRAPOLATE). If the preference
% CHEBTECH.PREF.EXTRAPOLATE is TRUE, then the values at the endpoints -1 and +1
% are always extrapolated (i.e., regardless of whether they evaluate to NaN).
%
% The CHEBTECH1 class supports the representation of vector-valued functions
% (for example, f = chebtech1(@(x) [sin(x), cos(x)])). In such cases, the values
% and coefficients are stored in a matrix (column-wise), and each
% component of the multi-valued function is truncated to the same length, even
% if the demands of 'happiness' imply that one of the components could be
% truncated to a shorter length than the others. All CHEBTECH1 methods should
% accept vectorised forms. Note that this representation is distinct from
% an array of CHEBTECH1 objects, such as [chebtech1(@(x) sin(x), chebtech1(@(x)
% cos(x)], for which there is little to no support.
%
% Class diagram: [<<CHEBTECH>>] <-- [CHEBTECH1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Methods implemented by this m-file:
    methods
        
        function obj = chebtech1(op, vscale, hscale, pref)
            % Constructor for the CHEBTECH1 class.
            
            % Return an empty CHEBTECH1 on null input:
            if ( (nargin == 0) || isempty(op) )
                return
            end
           
            % Define vscale if none given:
            if ( (nargin < 2) || isempty(vscale) )
                vscale = 0;
            end
            % Define hscale if none given:
            if ( (nargin < 3) || isempty(hscale) )
                hscale = 1;
            end
            % Determine preferences if not given, merge if some are given:
            if ( (nargin < 4) || isempty(pref) )
                pref = chebtech.pref;
            else
                pref = chebtech.pref(pref);
            end
            
            % Force nonadaptive construction if PREF.CHEBTECH.N is numeric:
            if ( ~isempty(pref.chebtech.n) && ~isnan(pref.chebtech.n) )
                % Evaluate op on the Chebyshev grid of given size:
                op = feval(op, chebtech1.chebpts(pref.chebtech.n));
            end
            
            % Actual construction takes place here:
            obj = populate(obj, op, vscale, hscale, pref);
            
            % In cases other than these, we will check for NaN values:
            if ( obj.ishappy || iscell(op) || isnumeric(op) )
                return
            end
            
            % Check for NaNs (if not happy):
            if ( any(isnan(obj.values(:))) )
                error('CHEBFUN:CHEBTECH1:constructor:naneval', ...
                    'Function returned NaN when evaluated.')
            end
            
        end
        
    end
    
    %% Static methods implemented by CHEBTECH1 class.
    % (This list is alphabetical.)
    methods ( Static = true )
        
        % Aliasing:
        coeffs = alias(coeffs, m)
        
        % Evaluate a Chebyshev interpolant using 2nd form barycentric formula:
        out = bary(x, values)
        
        % Compute Chebyshev barycentric weights:
        w = barywts(n)

        % Convert values to coefficients:
        coeffs = chebpoly(values)
        
        % Convert coefficients to values:
        values = chebpolyval(coeffs)
        
        % Compute Chebyshev points (x) and optionally quadrature (w)
        % and barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % Make a CHEBTECH1 (constructor shortcut):
        f = make(varargin);
        
        % Refinement function for CHEBTECH1 construction (evaluates OP on grid):
        [values, opints, giveUp] = refine(op, values, pref)
        
        % Compute Chebyshev quadrature weights:
        w = quadwts(n)
        
    end
    
end
