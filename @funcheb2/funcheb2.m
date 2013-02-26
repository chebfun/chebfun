classdef funcheb2 < funcheb
%FUNCHEB2   Approximate smooth functions on [-1,1] with Chebyshev interpolants. 
%
%   Class for approximating smooth functions on the interval [-1,1] using
%   function values at 2nd-kind Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   FUNCHEB2(OP) constructs a FUNCHEB2 object from the function handle OP. OP
%   should be vectorised (i.e., accept a vector input) and ouput a vector of the
%   same length. FUNCHEB2 objects allow for vectorised construction (i.e., of a
%   multi-valued function), in which case OP should accept a column vector of
%   length N and return a matrix of size NxM.
%
%   FUNCHEB2(OP, VSCALE) constructs a FUNCHEB2 with 'happiness' (see
%   HAPPINESSCHECK.m) relative to the maximum of the given vertical scale VSCALE
%   and the (column-wsie) infinity norm of the sampled function values of OP. If
%   not given, the VSCALE defaults to 0 initially.
%
%   FUNCHEB2(OP, VSCALE, HSCALE) uses a 'happiness' to both the vertical scale
%   VSCALE (as above) and the horizontal scale HSCALE. If not given, this
%   defaults to 1.
%
%   FUNCHEB2(OP, VSCALE, HSCALE, PREF) overrides the default behavior with that
%   given by the preference structure PREF. The constructor will also accept
%   inputs of the form FUNCHEB2(OP, PREF), but this usage is not advised.
%   Similarly, one can pass FUNCHEB2(OP, VSCALE, PREF). Furthermore, one can
%   replace PREF by TOL, the desired tolerance of the construction, which is
%   equivelent to passing a PREF with PREF.FUNCHEB2.eps = TOL.
%
%   FUNCHEB2(VALUES, ...) returns a FUNCHEB2 object which interpolates the
%   values in the columns of VALUES at 2nd-kind Chebyshev points and
%   FUNCHEB2({VALUES, COEFFS}, ... ) uses the Chebyshev coefficients passed in
%   COEFFS rather than computing them.
%
% Examples:
%   % Basic construction:
%   f = funcheb2(@(x) sin(x))
%
%   % Construction with preferences:
%   p = funcheb2.pref('sampletest', 0); % See help('funcheb2.pref') for details
%   f = funcheb2(@(x) sin(x), p)
%
%   % Vector-valued construction:
%   f = funcheb2(@(x) [sin(x), cos(x), exp(x)])
%
% See also FUNCHEB2.pref, FUNCHEB2.chebpts, FUNCHEB2.happinesscheck, FUNCHEB2.refine.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCHEB2 Class Description:
%
% The FUNCHEB2 class represents smooth functions on the interval [-1,1] using
% function values at 2nd-kind Chebyshev points and coefficients of the
% corresponding first-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 2nd-kind grid (see REFINE.m) until the
% representation is deemed 'happy' (see HAPPINESSCHECK.m). The resulting object
% can be used to evaluate and operate on the input function.
%
% The vertical scale VSCALE is used to enforce scale invariance in FUNCHEB2
% construction and subsequent operations. For example, that 
%          funcheb2(@(x) 2^300*f(x)) = 2^300*funcheb2(@(x) f(x)). 
% VSCALE may be optionally passed during to the constructor (if not, it defaults
% to 0), and during construction it is updated to be the maximum magnitude of
% the sampled function values.
%
% If the input operator OP evaluates to NaN or Inf at any of the sample points
% used by the constructor, then a suitable replacement is found by extrapolating
% (globally) from the numeric values (see EXTRAPOLATE.M). If the preference
% funcheb2.pref('extrapolate', TRUE) is set, then the endpoint values -1 and +1
% are always extrapolated (i.e., regardless of whether they evaluate to NaN).
%
% The FUNCHEB2 class supports the representation of vector-valued functions
% (for example, f = funcheb2(@(x) [sin(x), cos(x)])). In such cases, the values
% and coefficients are stored in a matrix (column-wise), and as such each
% component of the multi-valued function is truncated to the same length, even
% if the demands of 'happiness' imply that one of the components could be
% truncated to a shorter length than the others. All FUNCHEB2 methods should
% accept such vectorised forms. Note that this representation is distinct from
% an array of funcheb2 objects, such as [funcheb2(@(x) sin(x), funcheb2(@(x)
% cos(x)], for which there is little to no support.
%
% Class diagram: [<<funcheb>>] <-- [FUNCHEB2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Methods implemented by this m-file.
    methods
        
        function obj = funcheb2(op, vscale, hscale, pref)
            %Constructor for the FUNCHEB2 class.
            
            % Return an empty FUNCHEB2 on null input:
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
                pref = funcheb2.pref(vscale);
                vscale = 0;
                hscale = 1;
            elseif ( nargin == 3 && isstruct(hscale) )
                % hscale was actually a preference.
                pref = funcheb2.pref(hscale);
                hscale = 1;                
            elseif ( nargin < 4 )
                % Create:
                pref = funcheb2.pref;
            elseif ( ~isstruct(pref) )
                % An eps was passed.
                pref = funcheb2.pref('eps', pref);
            else
                % Merge:
                pref = funcheb2.pref(pref);
            end

            % Force nonadaptive construction if pref.funcheb2.n is numeric:
            if ( ~isempty(pref.funcheb2.n) && ~isnan(pref.funcheb2.n) )
                % Evaluate the op on the Chebyshev grid of given size:
                op = feval(op, funcheb2.chebpts(pref.funcheb2.n));
            end
            
            % Actual construction takes place here:
            obj = populate(obj, op, vscale, hscale, pref);
            
            % Check for NaNs: (if not happy)
            if ( ~obj.ishappy && pref.funcheb2.extrapolate )
                % Check for NaNs in interior only: (because extrapolate was on!)
                if ( any(any(isnan(obj.values(2:end-1,:)))))
                    error('CHEBFUN:FUNCHEB2:constructor:naneval', ...
                        'Function returned NaN when evaluated.')
                end
                % We make sure not to return NaNs at +1 and -1.
                valuesTemp = funcheb2.extrapolate(obj.values);
                obj.values([1, end],:) = valuesTemp([1, end],:);
            elseif ( ~obj.ishappy )
                % Here we throw an error if NaNs were encountered anywhere.
                if ( any(isnan(obj.values(:))) )
                    error('CHEBFUN:FUNCHEB2:constructor:naneval2', ...
                        'Function returned NaN when evaluated.')
                end
            end
            
        end
        
    end
    
    %% Static methods implemented by FUNCHEB2 class.
    % (This list is alphabetical)
    methods ( Static = true )
        
        % Evaluate a Chebyshev interpolant using barycentric formula.
        out = bary(x, values, kind)
        
        % Compute Chebyshev barycentric weights.
        w = barywts(n)
        
        % Convert values to coefficients.
        coeffs = chebpoly(values)
        
        % Convert coefficients to values.
        values = chebpolyval(coeffs)
        
        % Compute Chebyshev points of first kind (x) and optionally quadrature
        % (w) and barycentric (v) weights.
        [x, w, v] = chebpts(n)
        
        % Extrapolate (for NaNs / endpoints).
        [values, maskNaN, maskInf] = extrapolate(values)
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Refinement function for FUNCHEB2 construction. (Evaluates OP on grid)
        [values, opints, giveUp] = refine(op, values, pref)
        
        % Compute Chebyshev quadrature weights.
        w = quadwts(n)
        
        % Test the FUNCHEB2 class.
        pass = test(varargin);
        
        % Make a FUNCHEB2. (Constructor shortcut)
        f = make(varargin);
    end
    
    %% Private static methods implemented by FUNCHEB2 class.
    methods ( Access = public, Static = true )
        
        %
        
    end
    
end
