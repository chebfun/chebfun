function nVars = numVars(N)
%NUMVARS   Detect how many unknown variables a CHEBOP operates on.
%   NVARS = NUMVARS(N), where N is a CHEBOP, returns the integer NVARS which is
%   equal to the number of variables that N operates on.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% NUMVARS indicate how many unknown function we seek. This can be tricky of the
% operator is specified on CHEBMATRIX syntax, e.g. via
%   N.op = @(x, u) [diff(u{1}) + u{2}; u{1} + diff(u{2})];
% is in this case, nargin(N.op) does not match the number of variables we need.
% If nargin(N.op) is greater than two, we can however safely assume that N.op is
% specified on the form
%   N.op = @(x, u, v) [diff(u) + v; u + diff(v)];

% We begin by looking at whether we have the easy case!
narginN = nargin(N);
if ( narginN > 2 )
    % Need to subtract 1 since x is the first argument:
    nVars = narginN - 1;
    
elseif ( ~isempty(N.numVars) )
    % Now we know we're dealing with
    %   N.op = @(x, u) ...
    % But we don't know yet whether this is a system or not. In an ideal world,
    % the user has passed this information through the numVars property of the
    % CHEBOP.
    
    nVars = N.numVars; % Lucky us!
    
else
    % If that field is empty, we try to inspect the string representation of
    % N.op, and look for the highest index appearing inside {} -- this will only
    % work if N.op is specified as an anonymous function, not if it is a handle
    % to another function.

    % Obtain the function string:
    NopString = func2str(N.op);

    % Try to find out what the unknown function is. But first, we must check
    % whether we actually have the argument list available to us...
    firstRightPar = min(strfind(NopString, ')'));
    if ( isempty(firstRightPar) )
        % Don't have a list of arguments available. Take numVars == 1, and hope
        % for the best...
        warning('CHEBFUN:CHEBOP:numVars:numberOfArguments', ...
            ['Unable to determine the number of variable that the ', ...
            'chebop\noperates on. Assuming problem is a scalar problem,',...
            ' results might be\nunreliable. Please specify the number', ...
            ' of variables that the operator\n', ...
            'operates on via CHEBOP.numVars.'])
        nVars = 1;

    else
        % If nargin(N) == 2, the first variable appearing will be the
        % independent space variable, but if nargin(N) == 1, we only have the
        % unknown function in the argument list.
        if ( narginN == 1 )
            % The unknown variable appears between @( and the first ).
            firstLeftPar = min(strfind(NopString, '('));
            variableName = NopString(firstLeftPar+1 : firstRightPar - 1);

        else
            % The unknown variable appears between the first , and the first ),
            % e.g. @(x, ___ ).
            firstLeftComma = min(strfind(NopString, ','));
            variableName = NopString(firstLeftComma+1 : firstRightPar - 1);
        end

        % The regular expression we seek must include the variable name:
        expression = [variableName, '{[1-9]+}'];

        % Obtain all matches:
        match = regexp(NopString, expression, 'match');

        % If we don't have any matches, problem can't have been specified using
        % CHEBMATRIX syntax. But nargin(N) == 2, and we did indeed have the
        % argument list to the anonymous function passed, so we must be dealing
        % with a scalar problem!
        if ( isempty(match) )
            nVars = 1;
        else
            % Throw away the variable name and the { }, e.g. convert 'u{1}' to
            % '1':
            match = strrep(match, [variableName, '{'], '');
            match = strrep(match, '}', '');

            % We are now left with cell-array of strings that only contain
            % numbers. So convert to doubles!
            indx = str2double(match);

            % The number of variables that the CHEBOP operates on is the
            % greatest index that appears:
            nVars = max(indx);
        end
        
    end

end
