classdef (InferiorClasses = {?bndfun, ?unbndfun}) deltafun < fun
%DELTAFUN   Class for distributions based on Dirac-deltas on arbitrary intervals
%
%   Class for approximating generalized functions on the interval [a, b].
%   The smooth or classical part of the function is approximated by a
%   CLASSICFUN object while the Dirac delta functions are represented by a
%   pair of properties: DELTAMAG and LOCATION, which store the magnitude
%   and the location of the delta functions respectively. DELTAMAG is
%   generally a matrix, with its first row representing the delta functions,
%   while derivatives of delta functions are represented by higher rows of
%   this matrix. LOCATION is a vector, each element of which corresponds to
%   the location of a column in the DELTAMAG matrix.
%
% Constructor inputs:
%   DELTAFUN([], DATA) creates a DELTAFUN object with an empty FUNPART, with
%   the delta functions and their locations specified by the following fields in
%   the DATA structure:
%     DATA.DELTAMAG    (Default:  Empty)
%         Delta function magnitude matrix.  (See description above.)
%     DATA.DELTALOC    (Default:  Empty)
%         Row vector of locations of delta functions.
%
%   DELTAFUN(FUNPART, DATA) creates a DELTAFUN object with FUNPART as its
%   smooth function, while the locations and magnitudes of the delta functions
%   are specified by the DATA structure as above.
%
%   DELTAFUN(FUNPART, DATA, PREF) is the same as above but uses PREF to pass
%   any preferences.
%
%   DELTAFUN(OP, DATA, PREF) will call the CLASSICFUN constructor to create a
%   FUNPART out of the operator OP.  Fields in DATA which are not among those
%   listed above as recognized by DELTAFUN will be passed to the CLASSICFUN
%   constructor as-is.
%
%   DELTAFUN(F, ...) where F is a DELTAFUN simply returns F, i.e., other inputs
%   are ignored.
%
% See also CLASSICFUN, ONEFUN, FUN

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        
        % Smooth part of the representation.
        funPart     % (classical function which is a CLASSICFUN object)
        
        % DELTAMAG is a matrix storing signed magnitude of the delta functions
        % contained in a DELTAFUN object. The first row corresponds to the
        % delta functions themselves while higher order rows represent
        % derivatives of delta functions.
        deltaMag % double matrix [M X N]
        
        % Location of the delta functions.
        deltaLoc % double vector [1 X N]     
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        function obj = deltafun(op, data, pref)
            
            % Parse inputs.
            if ( (nargin <= 0) || isempty(op) )
                obj.funPart = [];
                obj.deltaMag = [];
                obj.deltaLoc = [];
                return
            end

            if ( (nargin < 2) || isempty(data) )
                data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = chebfunpref();
            else
                pref = chebfunpref(pref);
            end

            data = parseDataInputs(data, pref);

            % Given a DELTAFUN, return it.
            if ( isa(op, 'deltafun') )
                obj = op;
                return
            end

            % Assemble the funPart.
            if ( ~isa(op, 'classicfun') )
                op = classicfun.constructor(op, data, pref);
            end
            obj.funPart = op;

            % If there is no delta function information, we're done.
            if ( isempty(data.deltaMag) || isempty(data.deltaLoc) )
                if ( xor(isempty(data.deltaMag), isempty(data.deltaLoc)) )
                    warning('CHEBFUN:DELTAFUN:inconsistentData', ...
                        'Inconsistent deltaLoc and deltaMag.')
                end

                obj.deltaMag = [];
                obj.deltaLoc = [];
                return
            end

            % Make sure location is a row vector.
            if ( ~isempty(data.deltaLoc) )
                if ( min(size(data.deltaLoc)) > 1 )
                    error('CHEBFUN:DELTAFUN:deltafun:dim', ...
                        'deltaLoc should be a vector.');
                end
                data.deltaLoc = data.deltaLoc(:).';
            end

            % Check sizes:
            if ( ~isempty(data.deltaMag) && ...
                    (size(data.deltaMag, 2) ~= length(data.deltaLoc)) )
                error('CHEBFUN:DELTAFUN:deltafun:dim', ...
                    ['Impulse matrix (deltaMag) should have the same number' ...
                     ' of columns as locations (deltaLoc).']);
            end

            % Locations of delta functions should be within the domain:
            if ( ~isempty(data.deltaLoc) && ~isempty(obj.funPart) )
                dom = obj.funPart.domain;
                if ( (max(data.deltaLoc) > dom(2)) || ...
                     (min(data.deltaLoc) < dom(1)) )
                    error('CHEBFUN:DELTAFUN:deltafun:domain', ...
                        'Location of a delta fun is outside the domain.');
                end
            end

            % All checks done, assign inputs to the current object:
            obj.deltaMag = data.deltaMag;
            obj.deltaLoc = data.deltaLoc;

            % Simplify to merge redundant delta functions:
            obj = simplifyDeltas(obj, pref);
            if ( ~isa(obj, 'deltafun') )
                obj = deltafun(obj);
            end
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        % True if the DELTAFUN object has no delta functions       
        out = anyDelta(f)

        % Compose a DELTAFUN with an affine operator or another DELTAFUN
        f = compose(f, op, g, data, pref)
        
        % Complex conjugate of a DELTAFUN.
        f = conj(f)
        
        % DELTAFUN objects are not transposable.
        f = ctranspose(f)
        
        % Indefinite integral of a DELTAFUN.
        [F, jumpEnd] = cumsum(f, k, dim, shift)
        
        % DELTAFUN display of useful data
        data = dispData(f)
                
        % Derivative of a DELTAFUN.
        f = diff(f, k)
                
        % Evaluate a DELTAFUN.
        y = feval(f, x, varargin)
        
        % Flip columns of an array-valued DELTAFUN object.
        f = fliplr(f)
        
        % Flip/reverse a DELTAFUN object.
        f = flipud(f)
        
        % Imaginary part of a DELTAFUN.
        f = imag(f)
        
        % Compute the inner product of two DELTAFUN objects.
        out = innerProduct(f, g)
                
        % Always returns true, since DELTAFUN manages delta functions.
        out = isdelta(f)

        % True for an empty DELTAFUN.
        out = isempty(f)
        
        % Test if DELTAFUN objects are equal.
        out = isequal(f, g)

        % Test if a DELTAFUN is bounded.
        out = isfinite(f)
        
        % Test if a DELTAFUN is unbounded.
        out = isinf(f)
        
        % Test if a DELTAFUN has any NaN values.
        out = isnan(f)
        
        function out = isPeriodicTech(f)
        %ISPERIODICTECH    Test if the smooth part of f is is constructed with a 
        %basis of periodic functions. 
        
            % Calls ISPERIODICTECH on the CLASSICFUN part.
            out = isPeriodicTech(f.funPart);
        end
        
        % True for real DELTAFUN.
        out = isreal(f)
                
        % Ture if the DELTAFUN object is zero
        out = iszero(f)
        
        % Length of a DELTAFUN.
        len = length(f)
        
        % Wrap a cell around a single DELTAFUN.
        g = mat2cell(f, varargin)

        % Global maximum of a DELTAFUN on [-1,1].
        [maxVal, maxPos] = max(f)
        
        % Global minimum of a DELTAFUN on [-1,1].
        [minVal, minPos] = min(f)
        
        % Global minimum and maximum of a DELTAFUN on [-1,1].
        [vals, pos] = minandmax(f)
        
        % Subtraction of two DELTAFUN objects.
        f = minus(f, g)
        
        % Left matrix divide for DELTAFUN objects.
        X = mldivide(A, B)
        
        % Right matrix divide for a DELTAFUN.
        X = mrdivide(B, A)
        
        % Multiplication of DELTAFUN objects.
        f = mtimes(f, c)
        
        % Basic linear plot for DELTAFUN objects.
        varargout = plot(f, varargin)
        
        % Data for plotting a DELTAFUN
        data = plotData(f, g, h);

        % Addition of two DELTAFUN objects.
        f = plus(f, g)       
        
        % Dividing two DELTAFUNs
        f = rdivide(f, g)
        
        % Real part of a DELTAFUN.
        f = real(f)
        
        % Restrict a DELTAFUN to a subinterval.
        f = restrict(f, s)
        
        % Roots of a DELTAFUN.
        out = roots(f, varargin)
        
        % Size of a DELTAFUN.
        [siz1, siz2] = size(f, varargin)
        
        % Same locations in two DELTAFUNS
        varargout = sameDeltaLocs(f, g);
        
        % Simplify a DELTAFUN
        f = simplify(f, pref)
        
        % Simplify Delta functions
        f = simplifyDeltas(f, pref);
        
        % Definite integral of a DELTAFUN.
        out = sum(f, dim)
        
        % DELTAFUN multiplication.
        f = times(f, g)
        
        % Transfer delta function at the right end point to the next:
        [f, g] = transferDeltas(f, g);
        
        % DELTAFUN objects are not transposable.
        f = transpose(f)
        
        % Unary minus of a DELTAFUN.
        f = uminus(f)
        
        % Unary plus of a DELTAFUN.
        f = uplus(f)                
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % smooth fun constructor
        s = constructFunPart( op, pref)
        
        % remove zero columns
        [A, v] = cleanColumns(A, v, pref);
        
        % Create map
        map = createMap(ends)
        
        % remove zero trailing rows
        A = cleanRows(A, pref);
        
        % Find intersection based on some tolerance
        [x, idxV, idxW] = numIntersect( V, W, tol)
        
        % Constructor shortcut
        f = make(varargin)
        
        % Merge columns of a matrix based on duplicate values in v.
        [A, v] = mergeColumns(A, v, pref)
        
        % Merge delta function matrix
        [D, w] = mergeDeltas(A, v, B, u);
        
        % Costruct a zero DELTAFUN
        s = zeroDeltaFun(domain)
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHODS IMPLEENTED IN THIS FILE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = parseDataInputs(data, pref)
%PARSEDATAINPUTS   Parse inputs from the DATA structure and assign defaults.

if ( ~isfield(data, 'deltaMag') || isempty(data.deltaMag) )
    data.deltaMag = [];
end

if ( ~isfield(data, 'deltaLoc') || isempty(data.deltaLoc) )
    data.deltaLoc = [];
end

end
