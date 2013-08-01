classdef chebfun
%CHEBFUN    CHEBFUN class for representing piecewise smooth functions on [a, b].

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEBFUN Class Description:
%
% The CHEBFUN class is for representations of piecewise smooth functions on the
% interval [a, b].
%
% The CHEBFUN class is the main interface that a user will be using. Therefore,
% documentation should written appropriately, and checking of input arguments
% may be a little more forgiving.
%
% Class diagram: [jacfun] <>-- [CHEBFUN] <>-- [<<fun>>]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Access = public)
        domain
        funs
        impulses = [];
        isTransposed = 0;
    end

    methods

        function f = chebfun(varargin)
            % The main CHEBFUN constructor!

            % Return an empty CHEBFUN:
            if ( nargin == 0 || isempty(varargin{1}) )
                return
            end

            % Parse inputs:
            [op, dom, pref] = inputParser(varargin{:});

            if ( iscell(op) && all(cellfun(@(x) isa(x, 'fun'), op)) )
                % Construct a CHEBFUN from an array of FUN objects:

                if ( nargin > 1 )
                    error('CHEBFUN:chebfun:nargin', ...
                        'Only one input when passing an array of funs.')
                end
                
                f.funs = op;
                % Collect the domains together:
                dom = cellfun(@(fun) get(fun, 'domain'), f.funs, ...
                    'uniformOutput', false);
                f.domain = unique([dom{:}]);
                % Update values at jumps (first row of impulses):
                f.impulses = chebfun.jumpVals(f.funs, f.domain);

            else
                % Construct from function_handle, numeric, or string input:

                % Convert string input to function_handle:
                if ( ischar(op) )
                    op = str2op(op);
                end
                
                % Call the main constructor:
                [f.funs, f.domain] = chebfun.constructor(op, dom, pref);

                % Update values at jumps (first row of imps):
                f.impulses = chebfun.jumpVals(f.funs, f.domain, op);

                % Remove unnecessary breaks (but not those that were given):
                [ignored, index] = setdiff(f.domain, dom);
                f = merge(f, index', pref);

            end

        end

    end

    % Static methods implimented by CHEBFUN class.
    methods (Static = true)

        % Retrieve and modify preferences for this class.
        prefs = pref(varargin);

        % Splitting constructor.
        [funs, domain, op] = constructor(op, domain, pref);

        % Edge detector.
        [edge, vscale] = detectEdge(op, domain, hscale, vscale, derHandle);

        % Determine values of chebfun at breakpoints.
        vals = jumpVals(funs, ends, op);

    end
    
    % Static methods implimented by CHEBFUN class.
    methods (Static = true, Access = private)

        % Parse the inputs to the CHEBFUN construtor.
        [op, domain, pref] = inputParser(op, domain, varargin);
        
        % Convert a string input to a function_handle.
        op = str2op(op);

    end

    % Methods implimented by CHEBFUN class.
    methods

        % Display a CHEBFUN object.
        display(f);

        % Accuracy estimate of a CHEBFUN object.
        out = epslevel(f);
        
        % Retrieve and modify preferences for this class.
        out = get(f, prop);
        
        % Horizontal scale of a CHEBFUN object.
        out = hscale(f);
        
        % Vertical scale of a CHEBFUN object.
        out = vscale(f);

        % Size of a CHEBFUN object.
        [s1, s2] = size(f, dim);   
        
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                (Private) Methods implented in this m-file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function op = str2op(op)
    % This is here as it's a clean function with no other variables hanging
    % around in the scope.
    depVar = symvar(op);
    if ( numel(depVar) ~= 1 )
        error('CHEBFUN:STR2OP:indepvars', ...
            'Incorrect number of independent variables in string input.');
    end
    op = eval(['@(' depVar{:} ')', op]);
end

function [op, domain, pref] = inputParser(op, domain, varargin)
        % Parse inputs.

        args = varargin;
        if ( nargin == 1 )
            % chebfun(op)
            pref = chebfun.pref();
            domain = pref.chebfun.domain;
        elseif ( isstruct(domain) )
            % chebfun(op, pref)
            pref = domain;
            domain = pref.chebfun.domain;
        elseif ( nargin < 3 )
            % chebfun(op, domain)
            pref = chebfun.pref();
        elseif ( isstruct(varargin{1}) )
            % chebfun(op, domain, pref)
            pref = chebfun.pref(args{1});
            args(1) = [];
        elseif ( ~isnumeric(domain) )
            % chebfun(op, prop1, val1, ...)
            pref = chebfun.pref();
            args = [domain, args];
            domain = pref.chebfun.domain;
        else
            % chebfun(op, domain, prop1, val1, ...)
            pref = chebfun.pref;
        end

        % Obtain additional preferences:
        while ( ~isempty(args) )
            % Update these preferences:
            pref = chebfun.pref(pref, args{1}, args{2});
            args(1:2) = [];
        end

end


