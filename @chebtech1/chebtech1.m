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
%   the same length. CHEBTECH1 objects allow for array-valued construction
%   (i.e., of array-valued function), in which case OP should accept a vector
%   of length N and return a matrix of size NxM, where M is number of columns
%   of the multi -valued function.
%
%   CHEBTECH1(OP, DATA) constructs a CHEBTECH2 using the additional data
%   supplied in the DATA structure.  Fields currently recognized are:
%     DATA.VSCALE    (Default:  0)
%     DATA.HSCALE    (Default:  1)
%         The constructor builds a CHEBTECH1 with 'happiness' (see
%         HAPPINESSCHECK.m) relative to the maximum of the given vertical scale
%         DATA.VSCALE and the (column-wise) infinity norm of the sampled
%         function values of OP, and the fixed horizontal scale DATA.HSCALE.
%   If any fields in DATA are empty or not supplied, or if DATA itself is empty
%   or not supplied, appropriate default values are set.
%
%   CHEBTECH1(OP, DATA, PREF) overrides the default behavior with that given by
%   the preference structure PREF.
%
%   CHEBTECH1(VALUES, ...) returns a CHEBTECH1 object which interpolates the
%   values in the columns of VALUES at 1st-kind Chebyshev points and
%   CHEBTECH1({VALUES, COEFFS}, ... ) uses the Chebyshev coefficients passed in
%   COEFFS rather than computing them. If COEFFS are passed, the resulting
%   CHEBTECH1 is always deemed 'happy'.
%
% Examples: % Basic construction: f = chebtech1(@(x) sin(x))
%
%   % Construction with preferences:
%   p.sampleTest = 0; % See CHEBTECH.TECHPREF for details
%   f = chebtech1(@(x) sin(x), [], [], p)
%
%   % Array-valued construction:
%   f = chebtech1(@(x) [sin(x), cos(x), exp(x)])
%
% See also CHEBTECH, CHEBTECH.TECHPREF, CHEBPTS, HAPPINESSCHECK, REFINE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBTECH1 Class Description:
%
% The CHEBTECH1 class represents smooth functions on the interval [-1,1] using
% function values at 1st-kind Chebyshev points and coefficients of the
% corresponding 1st-kind Chebyshev series expansion.
%
% The constructor is supplied with a handle that evaluates a given function on
% an increasingly fine Chebyshev 1st-kind grid (see REFINE.m) until the
% representation is deemed 'happy' (see HAPPINESSCHECK). The resulting
% object can be used to evaluate and operate on the input function.
%
% More information can be found in the CHEBTECH class definition.
%
% Class diagram: [<<CHEBTECH>>] <-- [CHEBTECH1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        function obj = chebtech1(op, data, pref)
            % Parse inputs.
            if ( (nargin == 0) || isempty(op) )
                % Return an empty CHEBTECH1 on null input:
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
                op = feval(op, chebtech1.chebpts(pref.fixedLength));
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
                % Use contour integrals ("turbo" mode) if requested.
                if ( pref.useTurbo )
                    obj = constructorTurbo(obj, op, pref);
                end

                % No need to error check if we are happy:
                return
            end

            % Check for NaNs (if not happy):
            if ( any(isnan(obj.coeffs(:))) )
                % Here we throw an error if NaNs were encountered anywhere.
                error('CHEBFUN:CHEBTECH1:chebtech1:nanEval', ...
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
        t = angles(n);
        
        % Evaluate a Chebyshev interpolant using proper barycentric formula:
        out = bary(x, values)
        
        % Compute Chebyshev barycentric weights:
        w = barywts(n)
        
        % Compute Chebyshev points (x) and optionally quadrature (w)
        % and barycentric (v) weights.
        [x, w, v, t] = chebpts(n)
        
        % Tensor product grid of Chebyshev points in 1D, 2D or 3D:
        out = tensorGrid(N, dom)
        
        % Convert coefficients to values:
        values = coeffs2vals(coeffs)
        
        % Make a CHEBTECH1 (constructor shortcut):
        f = make(varargin);

        % Compute Chebyshev quadrature weights:
        w = quadwts(n)
        
        % Refinement function for CHEBTECH1 construction (evaluates OP on grid):
        [values, points, giveUp] = refine(op, values, pref)
        
        % Return the value-based discretization class which uses CHEBTECH1: 
        function disc = returnValsDisc()
            disc = @chebcolloc1;
        end
        
        % Convert values to coefficients:
        coeffs = vals2coeffs(values)
        
    end
        
end
