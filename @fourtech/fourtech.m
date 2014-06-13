classdef fourtech < smoothfun

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %% Properties of FOURTECH objects.
    properties ( Access = public )

        % Values of FOURTECH at equally spaced points from [-1,1).
        % For array-valued FOURTECH objects, each column
        % represents the interpolated values of a single function.
        values % (nxm double)

        % Coefficients are represented for the complex exponential form of 
        % the interpolant. The coefficients are stored in descending 
        % order so that c_{(N-1)/2} is
        % the first entry and c_{-(N-1)/2} is the last. For array-valued FOURTECH
        % objects, each column represents the coefficients of a single function.
        coeffs % (nxm double)

        % Vertical scale of the FOURTECH. This is a row vector storing the
        % magnitude of the largest entry in each column of VALUES. It is
        % convenient to store this as a property.
        vscale = 0 % (1xm double >= 0)

        % Horizontal scale of the FOURTECH. Although FOURTECH objects have in
        % principle no notion of horizontal scale invariance (since they always
        % live on [-1,1)), the input OP may have been implicitly mapped.
        % HSCALE is then used to enforce horizontal scale invariance in
        % construction and other subsequent operations that require it. It
        % defaults to 1 and is never updated.
        hscale = 1 % (scalar > 0)

        % Boolean value designating whether the FOURTECH is 'happy' or not. See
        % HAPPINESSCHECK.m for full documentation.
        ishappy % (logical)

        % Happiness level to which the FOURTECH was constructed (See
        % HAPPINESSCHECK.m for full documentation) or a rough accuracy estimate
        % of subsequent operations (See FOURTECH class documentation for
        % details).
        epslevel % (double >= 0)
        
        % Boolean value designating whether the FOURTECH represents a
        % real-valued function. This allows us to always return a real
        % result for things like evaluating a fourtech.
        isReal % (logical)
        
    end

    %% METHODS IMPLEMENTED BY THIS M-FILE:
    methods
        
        function obj = fourtech(op, data, pref)
            %Constructor for the FOURTECH class.

            % Parse inputs.
            if ( (nargin == 0) || isempty(op) )
                % Return an empty FOURTECH on null input:
                return
            end

            if ( (nargin < 2) || isempty(data) )
                    data = struct();
            end

            % [TODO]:  Preferences.
            if ( (nargin < 3) || isempty(pref) )
                pref = fourtech.techPref();
            else
                pref = fourtech.techPref(pref);
            end

            data = parseDataInputs(data, pref);

            % Force nonadaptive construction if PREF.NUMPOINTS is numeric:
            if ( ~isempty(pref.numPoints) && ~isnan(pref.numPoints) )
                % Evaluate op on the Fourier grid of given size:
                vals = feval(op, fourtech.fourpts(pref.numPoints));
                vals(1,:) = 0.5*(vals(1,:) + feval(op, 1));
                op = vals;
            end

            % Actual construction takes place here:
            obj = populate(obj, op, data.vscale, data.hscale, pref);
            
        end
        
    end
    
    %% STATIC METHODS IMPLEMENTED BY THIS CLASS:
    methods ( Static = true )
        
        % Aliasing:
        coeffs = alias(coeffs, m)
        
        % Evaluate a Fourier interpolant using the barycentric formula:
        out = bary(x, values)
        
        % Compute Fourier barycentric weights:
        w = barywts(n)
        
        % Compute Fourier points (x) and optionally quadrature (w)
        % and barycentric (v) weights:
        [x, w, v] = fourpts(n);
        
        % Convert coefficients to values:
        values = coeffs2vals(coeffs);
        
        % Make a FOURTECH (constructor shortcut):
        f = make( varargin );
        
        % Compute Fourier quadrature weights (trapezoidal rule):
        w = quadwts(n)
        
        % Refinement function for FOURTECH construction (evaluates OP on grid):
        [values, points, giveUp] = refine(op, values, pref)
        
        % Convert values to coefficients:
        coeffs = vals2coeffs(values)
                
        % Retrieve and modify preferences for this class.
        p = techPref(q)
        
    end
    
    %% METHODS IMPLEMENTED BY THIS CLASS:
    methods
        
        % Compose two FOURTECH objects or a FOURTECH with a function handle:
        h = compose(f, op, g, data, pref)
        
        % Plot (semilogy) the Chebyshev coefficients of a FOURTECH object.
        h = coeffsplot(f, varargin)

        % Get method:
        val = get(f, prop);
        
    end
    
    methods       
        
        function c = legpoly(f)
            error('CHEBFUN:FOURTECH:legpoly:notAvailable',...
                'Cannot convert a Fourier based chebfun to a Legendre Series. Try first converting f to a Chebyshev based chebfun');
        end
        
        function [f, rootsLeft, rootsRight] = extractBoundaryRoots(f, numRoots)
            error('CHEBFUN:FOURTECH:extractBoundaryRoots:notAvailable',...
                'Function not implemented. Try first converting f to a Chebyshev based chebfun');
        end
        
    end
    
end

function data = parseDataInputs(data, pref)
%PARSEDATAINPUTS   Parse inputs from the DATA structure and assign defaults.

if ( ~isfield(data, 'vscale') || isempty(data.vscale) )
    data.vscale = 0;
end

if ( ~isfield(data, 'hscale') || isempty(data.hscale) )
    data.hscale = 1;
end

end