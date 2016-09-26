classdef trigtech < smoothfun % (Abstract)
%TRIGTECH   Approximate smooth periodic functions on [-1,1] with trigonometric 
%           interpolants.
%
%   Class for approximating smooth periodic functions on the interval [-1,1]
%   using function values at equally spaced points on [-1,1).
%
% Constructor inputs:
%   TRIGTECH(OP) constructs a TRIGTECH object from the function handle OP. OP
%   should be vectorized (i.e., accept a vector input) and ouput a vector of the
%   same length. TRIGTECH objects allow for array-valued construction (i.e., of
%   array-valued function), in which case OP should accept a vector of length N
%   and return a matrix of size NxM, where M is number of columns of the multi
%   -valued function.
%
%   TRIGTECH(OP, DATA) constructs a TRIGTECH using the additional data
%   supplied in the DATA structure.  Fields currently recognized are:
%     DATA.VSCALE    (Default:  0)
%     DATA.HSCALE    (Default:  1)
%         The constructor builds a TRIGTECH with 'happiness' (see
%         HAPPINESSCHECK.m) relative to the maximum of the given vertical scale
%         DATA.VSCALE and the (column-wise) infinity norm of the sampled
%         function values of OP, and the fixed horizontal scale DATA.HSCALE.
%   If any fields in DATA are empty or not supplied, or if DATA itself is empty
%   or not supplied, appropriate default values are set.
%
%   TRIGTECH(OP, DATA, PREF) overrides the default behavior with that given by
%   the preference structure PREF.
%
%   TRIGTECH(VALUES, ...) returns a TRIGTECH object which interpolates the
%   values in the columns of VALUES at equally spaced points and
%   TRIGTECH({VALUES, COEFFS}, ... ) uses the trig-coefficients passed in
%   COEFFS rather than computing them. If COEFFS are passed, the resulting
%   TRIGTECH is always deemed 'happy'.
%
% Examples: % Basic construction: f = trigtech(@(x) exp(sin(pi*x)))
%
%   % Construction with preferences:
%   p.sampleTest = 0; % See TRIGTECH.TECHPREF for details
%   f = trigtech(@(x) sin(x), [], [], p)
%
%   % Array-valued construction:
%   f = trigtech(@(x) tanh([sin(pi*x), cos(pi*x), cos(pi*sin(pi*x))]))
%
% See also TRIGTECH.TECHPREF, TRIGPTS, HAPPINESSCHECK, REFINE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIGTECH Class Description:
%
% The TRIGTECH class is for representations of smooth periodic functions on the
% interval [-1,1] via interpolated function values at equally spaced points
% using trigonometric Fourier series.
%
% The vertical scale VSCALE is used to enforce scale invariance in TRIGTECH
% construction and subsequent operations. For example, that
%   trigtech(@(x) 2^300*f(x)) = 2^300*trigtech(@(x) f(x)).
%
% VSCALE may be optionally passed to the constructor (if not, it defaults to 0),
% and during construction it is updated to be the maximum magnitude of the
% sampled function values. Similarly the horizontal scale HSCALE is used to
% enforce scale invariance when the input OP has been implicitly mapped from a
% domain other than [-1 1] before being passed to the TRIGTECH constructor.
%
% If the input operator OP in a call to TRIGTECH evaluates to NaN or Inf at
% any of the sample points used by the constructor, then an error is thrown.
%
% The TRIGTECH class supports the representation of array-valued functions (for
% example, f = trigtech(@(x) [sin(pi*x), cos(pi*x)])). In such cases, the values
% and coefficients are stored in a matrix (column-wise), and as such each
% component of the array-valued function is truncated to the same length, even
% if the demands of 'happiness' imply that one of the components could be
% truncated to a shorter length than the others. All TRIGTECH methods should
% accept such array-valued forms. Note that this representation is distinct from
% an array of TRIGTECH objects, for which there is little to no support.
%
% Class diagram: [<<SMOOTHFUN>>] <-- [TRIGTECH]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )

        % Values of TRIGTECH at equally spaced points from [-1,1). For
        % array-valued TRIGTECH objects, each column represents the interpolated
        % values of a single function.
        values % (nxm double)

        % Coefficients are represented for the complex exponential form of the
        % interpolant. The coefficients are stored in descending order so that
        % c_{(N-1)/2} is the first entry and c_{-(N-1)/2} is the last. For
        % array-valued TRIGTECH objects, each column represents the coefficients
        % of a single function.
        coeffs % (nxm double)

        % Boolean value designating whether the TRIGTECH is 'happy' or not.
        % See HAPPINESSCHECK.m for full documentation.
        ishappy % (logical)

        % Boolean value designating whether the TRIGTECH represents a
        % real-valued function. This allows us to always return a real result
        % for things like evaluating a TRIGTECH.
        isReal % (logical)
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        function obj = trigtech(op, data, pref)
            %Constructor for the TRIGTECH class.

            % Parse inputs.
            if ( (nargin == 0) || isempty(op) )
                % Return an empty TRIGTECH on null input:
                return
            end

            if ( (nargin < 2) || isempty(data) )
                data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = trigtech.techPref();
            else
                pref = trigtech.techPref(pref);
            end

            data = trigtech.parseDataInputs(data, pref);

            % Force nonadaptive construction if PREF.FIXEDLENGTH is numeric:
            if ( ~(isnumeric(op) || iscell(op)) && ...
                    ~isempty(pref.fixedLength) && ~isnan(pref.fixedLength) )
                % Evaluate op on the equi-spaced grid of given size:
                vals = feval(op, trigtech.trigpts(pref.fixedLength));
                vals(1,:) = 0.5*(vals(1,:) + feval(op, 1));
                op = vals;
            end

            % Actual construction takes place here:
            obj = populate(obj, op, data, pref);
            
            % Set length of obj to PREF.FIXEDLENGTH (if it is non-trivial).
            if ( (isnumeric(op) || iscell(op)) && ...
                    ~isempty(pref.fixedLength) && ~isnan(pref.fixedLength) )
                obj = prolong(obj, pref.fixedLength);
            end
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        function out = isPeriodicTech(f)
        %ISPERIODICTECH   True for TRIGTECH.
            out = 1;
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Static = true )
        
        % Aliasing:
        coeffs = alias(coeffs, m)

        % Differentiation matrix in Fourier basis.
        D = diffmat(n, p)
        
        % Compute trigonometric points (x) and optionally quadrature (w)
        % and barycentric (v) weights:
        [x, w] = trigpts(n);
        
        % Tensor product grid of equispaced points in 1D, 2D or 3D:
        out = tensorGrid(N, dom)
        
        % Convert coefficients to values:
        values = coeffs2vals(coeffs);
        
        % Horner scheme for evaluation of a TRIGTECH
        y = horner(x, c, allReal)

        % Make a TRIGTECH (constructor shortcut):
        f = make(varargin);
        
        % Compute trigonometric quadrature weights (trapezoidal rule):
        w = quadwts(n)
        
        % Refinement function for TRIGTECH construction (evaluates OP on grid):
        [values, points, giveUp] = refine(op, values, pref)
        
         % Return the value-based discretization class which uses TRIGTECH: 
        function disc = returnValsDisc()
            disc = @trigcolloc;
        end
        
        % Retrieve and modify preferences for this class.
        p = techPref(q)

        % Convert values to coefficients:
        coeffs = vals2coeffs(values)

        % Parse inputs passed to the constructor via the DATA argument.
        data = parseDataInputs(data, pref)
    end
    
end
