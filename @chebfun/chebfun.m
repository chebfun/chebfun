classdef chebfun
%CHEBFUN  CHEBFUN class for representing piecewise functions on [a, b].

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
    end
    
    methods
        
        function f = chebfun(op, domain, varargin)
            % The main CHEBFUN constructor!
            
            % Return an empty CHEBFUN:
            if ( nargin == 0 || isempty(op) )
                return
            end 

            % Parse inputs:
            args = varargin;
            if ( nargin == 1 )
                % chebfun(op)
                pref = chebfun.pref;
                domain = pref.chebfun.domain;
            elseif ( isstruct(domain) )
                % chebfun(op, pref)
                pref = chebfun.pref(domain);
                domain = pref.chebfun.domain;
            elseif ( nargin < 3 )
                % chebfun(op, domain)
                pref = chebfun.pref;
            elseif ( isstruct(varargin{1}) )
                % chebfun(op, domain, pref)
                pref = chebfun.pref(args{1});
                args(1) = [];   
            elseif ( ~isnumeric(domain) )
                % chebfun(op, prop1, val1, ...)
                pref = chebfun.pref;    
                args = [domain  args];
                domain = pref.chebfun.domain;
            else
                % chebfun(op, domain, prop1, val1, ...)
                pref = chebfun.pref;
            end
            
            stats = false;
            % Obtain additional preferences:
            while ( ~isempty(args) )
                if ( strcmpi(args{1}, 'stats') )
                    stats = strcmpi(args{2}, 'on') || args{2} == 1;
                    args(1:2) = [];
                    continue
                end
                % Update these preferences:
                pref = chebfun.pref(pref, args{1}, args{2});
                args(1:2) = [];
            end

            % Construct a chebfun from an array of funs:
            if ( iscell(op) && all(cellfun(@(x) isa(x, 'fun'), op)) )
                if ( nargin > 1 )
                    error('CHEBFUN:chebfun:nargin', ...
                        'Only one input when passing an array of funs.')
                end
                f.funs = op;
                if ( ~isempty(op) )
                    f.domain = op{1}.domain;
                end
                return
            end

            % Unwrap should also convert strings to function_handles, but
            % for now we do it here.
            if ( ischar(op) )
                op = str2op(op); 
            end
            
            if ( stats )
                op = @(x) statop(x, op);
                op(NaN);
            end
            
            % Call the main constructor:
            [f.funs, f.domain] = chebfun.constructor(op, domain, pref);
            
            % Update values at jumps, first row of imps:
            f.impulses = chebfun.jumpvals(f.funs, f.domain, op);   
            
            if ( stats )
                op('foo');
            end
            
        end
        
    end
    
    % Static methods implimented by CHEBFUN class.
    methods (Static = true)
        
        % Retrieve and modify preferences for this class.
        prefs = pref(varargin)
        
        % Splitting constructor.
        [funs, domain] = constructor(op, domain, pref)
        
        % Edge detector.
        [edge, vscale] = detectedge(op, domain, hscale, vscale, pref, d)
        
        % Determine values of chebfun at breakpoints.
        vals = jumpvals(funs, ends, op)
        
    end
    
end

function op = str2op(op)
    % This is here as it's a clean function with no other variables hanging
    % around in the scope.
    depvar = symvar(op);
    if ( numel(depvar) ~= 1 )
        error('CHEBFUN:STR2OP:indepvars', ...
            'Incorrect number of independent variables in string input.');
    end
    op = eval(['@(' depvar{:} ')' op]);
end

function y = statop(x, op)
    % Allows statistics to be returns about construction.
    
    % Can only do this for a function handle.
    if ( ~isa(op, 'function_handle') )
        warning('CHEBFUN:STSTOP:notafunc', 'Unable to return statistics.')
        y = op(x);
        return
    end
    
    % Persistent variables:
    persistent store_x 
    persistent caller

    % initialise them to empty cells:
    if ( isempty(store_x) )
        store_x = {};
        caller = {};
    end

    if ( isnan(x) )
        % Passing a NaN resets:
        store_x = {};
        caller = {};
        
    elseif ( ischar(x) )
        % Passing a string displays data:
        numEvals = sum(cellfun(@numel, store_x));
        fprintf('Total number of function evaluations = %d\n', numEvals);
        x = cellfun(@numel, store_x,'uniformoutput', false);
        disp([caller.' x.'])
        x = sort(vertcat(store_x{:}));
        hist(x, 9), shg
        store_x = {};
        caller = {};
    else
        if ( isempty(x) )
            return
        end
        % Else the operator is evaluated:
        y = op(x);
        %  and inputs are stored:
        store_x = [store_x x];
        tmp = dbstack;
        tmp(1:2) = [];
        if ( strncmpi(tmp(1).name, '@(x)', 4) )
            tmp(1) = [];
        end
        caller = [caller tmp(1).name];
        
    end

end


