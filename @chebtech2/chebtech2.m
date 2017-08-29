classdef chebtech2 < chebtech
%CHEBTECH2   Approximate smooth functions on [-1,1] with Chebyshev interpolants.
%
%   Class for approximating smooth functions on the interval [-1,1] using
%   function values at 2nd-kind Chebyshev points and coefficients of the
%   corresponding 1st-kind Chebyshev series expansion.
%
% Constructor inputs:
%   CHEBTECH2(OP) constructs a CHEBTECH2 object from the function handle OP. OP
%   should be vectorised (i.e., accept a vector input) and ouput a vector of
%   the same length. CHEBTECH2 objects allow for array-valued construction
%   (i.e., of an array-valued function), in which case OP should accept a column
%   vector of length N and return a matrix of size NxM.
%
%   CHEBTECH2(OP, DATA) constructs a CHEBTECH2 using the additional data
%   supplied in the DATA structure.  Fields currently recognized are:
%     DATA.VSCALE    (Default:  0)
%     DATA.HSCALE    (Default:  1)
%         The constructor builds a CHEBTECH2 with 'happiness' (see
%         HAPPINESSCHECK.m) relative to the maximum of the given vertical scale
%         DATA.VSCALE and the (column-wise) infinity norm of the sampled
%         function values of OP, and the fixed horizontal scale DATA.HSCALE. If
%         not given (or given as empty), the VSCALE defaults to 0 initially,
%         and HSCALE defaults to 1.
%   If any fields in DATA are empty or not supplied, or if DATA itself is empty
%   or not supplied, appropriate default values are set.
%
%   CHEBTECH2(OP, DATA, PREF) overrides the default behavior with that given by
%   the preference structure PREF. See CHEBTECH.TECHPREF for details.
%
%   CHEBTECH2(VALUES, ...) returns a CHEBTECH2 object which interpolates the
%   values in the columns of VALUES at 2nd-kind Chebyshev points and
%   CHEBTECH2({VALUES, COEFFS}, ... ) uses the Chebyshev coefficients passed in
%   COEFFS rather than computing them from VALUES. If VALUES, is empty then the
%   CHEBTECH2 is constructed directly from the COEFFS. If COEFFS are passed,
%   the resulting CHEBTECH2 is always deemed 'happy'.
%
% Examples: % Basic construction: f = chebtech2(@(x) sin(x))
%
%   % Construction with preferences:
%   p.sampleTest = 0; % See CHEBTECH.TECHPREF for details
%   f = chebtech2(@(x) sin(x), [], [], p)
%
%   % Array-valued construction:
%   f = chebtech2(@(x) [sin(x), cos(x), exp(x)])
%
% See also CHEBTECH, CHEBTECH.TECHPREF, CHEBPTS, HAPPINESSCHECK, REFINE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBTECH2 Class Description:
%
% The CHEBTECH2 class represents smooth functions on the interval [-1,1] using
% function values at 2nd-kind Chebyshev points and coefficients of the
% corresponding 1st-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 2nd-kind grid (see REFINE.m) until the
% representation is deemed 'happy' (see HAPPINESSCHECK.m). The resulting object
% can be used to evaluate and operate on the input function.
%
% More information can be found in the CHEBTECH class definition.
%
% Class diagram: [<<CHEBTECH>>] <-- [CHEBTECH2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function obj = chebtech2(op, data, pref)
            % Parse inputs.
            if ( (nargin == 0) || isempty(op) )
                % Return an empty CHEBTECH2 on null input:
                return
            end

            if ( (nargin < 2) || isempty(data) )
                    data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = chebtech.techPref();
            else
                pref = chebtech.techPref(pref);
            end

            data = chebtech.parseDataInputs(data, pref);

            % Force nonadaptive construction if PREF.FIXEDLENGTH is numeric and
            % we're not using contour integrals.
            if ( ~(isnumeric(op) || iscell(op)) && ...
                    ~isnan(pref.fixedLength) && ~pref.useTurbo )
                % Evaluate op on the Chebyshev grid of given size:
                op = feval(op, chebtech2.chebpts(pref.fixedLength));
            end

            % Actual construction takes place here:
            [obj, values] = populate(obj, op, data, pref);

            if ( isnumeric(op) || iscell(op) )
                % Set length of obj to PREF.FIXEDLENGTH (if it is non-trivial).
                if ( ~isnan(pref.fixedLength ) )
                    obj = prolong(obj, pref.fixedLength);
                end

                % No need to error check when constructing from discrete data.
                return
            elseif ( obj.ishappy )
                % Use contour integrals if requested.
                if ( pref.useTurbo )
                    obj = constructorTurbo(obj, op, pref);
                end

                % No need to error check if we are happy:
                return
            end

            % Check for NaNs (if not happy):
            if ( pref.extrapolate )
                % Check for NaNs in interior only (because extrapolate was on):
                if ( any(any(isnan(obj.coeffs(2:end-1,:)))) )
                    error('CHEBFUN:CHEBTECH2:chebtech2:nanEval', ...
                        'Function returned NaN when evaluated.')
                end
                % We make sure not to return NaNs at +1 and -1.
                valuesTemp = extrapolate(obj, values);
                valuesTemp([1 end], :) = valuesTemp([1,end],:);
                obj.coeffs = obj.vals2coeffs(valuesTemp);
            elseif ( any(isnan(obj.coeffs(:))) )
                % Here we throw an error if NaNs were encountered anywhere.
                error('CHEBFUN:CHEBTECH2:chebtech2:nanEval2', ...
                    'Function returned NaN when evaluated.')
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Aliasing:
        coeffs = alias(coeffs, m)
        
        % Angles of Chebyshev points. (i.e., acos(chebpts(n))
        t = angles(n)
        
        % Evaluate a Chebyshev interpolant using the barycentric formula:
        out = bary(x, values)
        
        % Compute Chebyshev barycentric weights:
        w = barywts(n)
        
        % Compute Chebyshev points (x) and optionally quadrature (w)
        % and barycentric (v) weights:
        [x, w, v, t] = chebpts(n);
        
        % Tensor product grid of Chebyshev points in 1D, 2D or 3D:
        out = tensorGrid(N, dom)
        
        % Convert coefficients to values:
        values = coeffs2vals(coeffs);
        
        % Make a CHEBTECH2 (constructor shortcut):
        f = make(varargin);
        
        % Compute Chebyshev quadrature weights:
        w = quadwts(n)
        
        % Refinement function for CHEBTECH2 construction (evaluates OP on grid):
        [values, points, giveUp] = refine(op, values, pref)
       
        % Return the value-based discretization class which uses CHEBTECH2: 
        function disc = returnValsDisc()
            disc = @chebcolloc2;
        end
        
        % Convert values to coefficients:
        coeffs = vals2coeffs(values)
        
    end
            
end
