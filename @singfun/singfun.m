classdef (InferiorClasses = {?chebtech2, ?chebtech1}) singfun < onefun %(See Notes)
%SINGFUN   Class for functions with singular endpoint behavior.
%
%   Class for approximating singular functions on the interval [-1,1] using a
%   smooth part with no singularities and two singular factors (1+x).^a and
%   (1-x).^b, where a and b are real numbers.
%
%   SINGFUN class description
%
%   The SINGFUN class represents a function of the form
%
%          f(x) = s(x) (1+x)^a (1-x)^b
%
%   on the interval [-1,1]. The exponents a and b are assumed to be real. The
%   constructor is supplied with a handle that evaluates the function f at any
%   given points within the interval [-1,1]. The endpoint values are likely to
%   return Inf or NaN results.
%
%   Ideally, the "smooth" function s(x) is analytic, or at least much more
%   compactly represented than f is. The resulting object can be used to
%   evaluate and operate on the function f. If a and b are unknown at the time
%   of construction, the constructor will try to determine appropriate values
%   automatically by sampling the function handle. Note, however, that this
%   process is not completely robust, and the singular terms in general do not
%   perfectly factor out singular behavior. The constructor can be forced to
%   consider only integer exponents.
%
%   Multiplication and division are as good as the corresponding operations on
%   the smooth part. Addition and subtraction are much less reliable, as the sum
%   of two SINGFUN objects with different exponents is not necessarily a
%   SINGFUN, nor a smooth function. If all but integer exponents can be factored
%   out of the summands, the process is fine, but in other circumstances the
%   process may throw an error.
%
%   SINGFUN(OP) constructs a SINGFUN object from the function handle OP. It
%   first tries to determine the type and the order of the singularities or the
%   roots at the endpoints by sampling near -1 and 1 and then forms a new
%   operator which governs the smooth part of OP by factoring out the end point
%   singularities. Finally it constructs an approximation of the smooth part by
%   calling the SMOOTHFUN constructor. The type and the order of the
%   singularities along with the SMOOTHFUN representation of the smooth part are
%   stored in corresponding member fields of the instantiation.
%
%   SINGFUN(OP, DATA) constructs a SINGFUN object using the data in the MATLAB
%   structure DATA.  DATA fields recognized by SINGFUN are.
%     DATA.SINGTYPE    (Default:  Determined by CHEBFUNPREF)
%         A 1x2 cell array of strings and type of the singularities specified
%         by these strings may help the singularity detector to determine the
%         order of the singularities more efficiently and save some computing
%         time.  Valid strings for SINGTYPE are 'none', 'pole', 'sing' or
%         'root'.
%     DATA.EXPONENTS   (Default:  Empty)
%         If DATA.EXPONENTS is nonempty, then instead of determining the
%         singularity types and orders by sampling the function values at -1
%         and 1, the constructor takes the values saved in the 1x2 vector
%         EXPONENTS as the orders of the singularities.
%   If any fields in DATA are empty or not supplied, or if DATA itself is empty
%   or not supplied, appropriate default values are set.  Any fields in DATA
%   which are not recognized will be passed as-is to the SMOOTHFUN constructor.
%
%   SINGFUN(OP, DATA, PREF) constructs a SINGFUN using the preferences given by
%   PREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEVELOPER NOTES:
%   It is a Matlab requirement to specify exactly which classes are inferior to
%   a given class. One might think that writing "InferiorClasses = {?smoothfun}"
%   should be OK but it turns out that subclasses do not inherit the attribute
%   of being inferior and all inferior clases/subclasses should be mentioned 
%   explicitly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % Smooth part of the representation.
        smoothPart      % (SMOOTHFUN)
        
        % Exponents of the singularities at the two endpoints.
        exponents       % (1x2 double)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function obj = singfun(op, data, pref)
            
            % Parse inputs.
            if ( nargin < 1 )
                % No input arguments; return an empty SINGFUN.
                obj.smoothPart = [];
                obj.exponents = [];
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

            % Leave SINGFUN OPs alone, and upgrade SMOOTHFUN OPs to SINGFUNs.
            if ( nargin == 1 )
                if ( isa(op, 'singfun') )
                    obj = op;
                    return
                elseif ( isa(op, 'smoothfun') )
                    obj.smoothPart = op;
                    obj.exponents = [0, 0];
                    return
                end
            end

            % Find exponents.poly
            if ( ~isempty(data.exponents) )
                % Check values of supplied exponents.
                if ( any(size(data.exponents) ~= [1, 2]) || ...
                        ~isa(data.exponents, 'double') )
                    error('CHEBFUN:SINGFUN:singfun:badExponents', ...
                        'Exponents must be a 1x2 vector of doubles.');
                end

                % If any of EXPONENTS is NaN, try to determine:
                maskNaN = isnan(data.exponents);
                if ( any(maskNaN) )
                    tmpExps = singfun.findSingExponents(op, data.singType);
                    data.exponents(maskNaN) = tmpExps(maskNaN);
                end

                obj.exponents = data.exponents;
            else
                % Exponents not given.  Determine them.
                obj.exponents = singfun.findSingExponents(op, data.singType);
            end

            % Make sure that op is a function handle, a smoothfun, or numeric.
            if ( ~isnumeric(op) && ~isa(op, 'function_handle') && ...
                    ~isa(op, 'smoothfun') )
                error( 'CHEBFUN:SINGFUN:singfun:badOp', ...
                    ['First argument must be a function handle or a ', ...
                     'SMOOTHFUN, not a %s.'], class(op));
            end

            % Check to avoid array-valued operators.
            if ( (isnumeric(op) && size(op, 2) > 1) || ...
                 (~isnumeric(op) && size(feval(op, 0), 2) > 1) )
                error('CHEBFUN:SINGFUN:singfun:arrayValued', ...
                    'SINGFUN does not support array-valued construction.');
            end

            % Extrapolate when the given function blows up.
            if ( any(obj.exponents < 0) )
                pref.techPrefs.extrapolate = true;
            end

            % Construct the smoothPart.
            if ( isa(op, 'smoothfun') )
                % smoothPart was handed to us.
                obj.smoothPart = op;
            elseif ( isnumeric(op) )
                % Constructing from numeric. No function handle to moidify.
                obj.smoothPart = singfun.constructSmoothPart(op, data, pref);
            else               
                % Loosen tolerance:
                if ( any(obj.exponents) )
                    pref.chebfuneps = max(pref.chebfuneps, 1e-14);
                end
                % Construct new function Handle by factorinf out singular terms
                % from the operator based on the values in exponents.
                smoothOp = singOp2SmoothOp(op, obj.exponents);
                obj.smoothPart = ...
                    singfun.constructSmoothPart(smoothOp, data, pref);
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )

        function out = isPeriodicTech(f)
        %ISPERIODICTECH   Test if the smooth part of f is is constructed with a
        %basis of periodic functions.
            out = isPeriodicTech(f.smoothPart);
        end

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % SmoothPart constructor
        s = constructSmoothPart( op, vscale, hscale, pref )
        
        % Method for finding the order of singularities
        exponents = findSingExponents( op, singType )
        
        % Find integer order singularities, i.e. poles
        poleOrder = findPoleOrder( op, SingEnd )
        
        % Finding fractional order singularities (+ve or -ve).
        branchOrder = findSingOrder( op, SingEnd )
        
        % Make a SINGFUN (constructor shortcut):
        f = make(varargin);
      
        % Convert a SMOOTHFUN to a SINGFUN.
        f = smoothFun2SingFun(f) 
        
        % Construct a zero SINGFUN
        s = zeroSingFun()
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OTHER FUNCTIONS IMPLEMENTED IN THIS M-FILE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function op = singOp2SmoothOp(op, exponents)
%SINGOP2SMOOTHOP   Converts a singular operator to a smooth one by removing the
%   singularity(ies).
%   SINGOP2SMOOTHOP(OP, EXPONENTS) returns a smooth operator OP by removing the
%   singularity(ies) at the endpoints -1 and 1. EXPONENTS are the order of the
%   singularities.
%
%   For example, opNew = singOp2SmoothOp(opOld, [-a -b]) means
%   opNew = opOld.*(1+x).^a.*(1-x).^b
%
% See also SINGFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( all(exponents) )
    % Both exponents are non trivial:
    op = @(x) op(x)./((1 + x).^(exponents(1)).*(1 - x).^(exponents(2)));
    
elseif ( exponents(1) )
    % (1+x) factor at the left end point:
    op = @(x) op(x)./(1 + x).^(exponents(1));
    
elseif ( exponents(2) )
    % (1-x) factor at the right end point:
    op = @(x) op(x)./(1 - x).^(exponents(2));
    
end

end

function data = parseDataInputs(data, pref)
%PARSEDATAINPUTS   Parse inputs from the DATA structure and assign defaults.

if ( ~isfield(data, 'exponents') || isempty(data.exponents) )
    data.exponents = [];
end

if ( ~isfield(data, 'singType') || isempty(data.singType) )
    defaultSingType = pref.blowupPrefs.defaultSingType;
    data.singType = {defaultSingType, defaultSingType};
end

end
