function guess = constructInit(initInput, allVarNames, indVarNameSpace, xt)
%CONSTRUCTINIT    Convert initial condition input from CHEBGUI to a CHEBFUN
%
%  Calling sequence:
%       GUESS = CONSTRUCTINIT(INITINPUT, ALLVARNAMES, INDVARNAMESPACE, XT)
%
%  where the inputs are:
%
%       INITINPUT:        The input in the initial guess/condition field of 
%           CHEBGUI.
%       ALLVARNAMES:      A cellstring containing the name of all variables that 
%           appear in the problem.
%       INDVARNAMESPACE:  The name of the independent space/time variable in the
%           problem.
%       XT:               The independent space/time variable on the domain that
%           the problem specifies.
%
%  and the output is:
%       GUESS: A CHEBFUN/CHEBMATRIX representing the initial guess/condition for
%           the problem.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% In order to be able to evaluate the input below, we need to assign r, x or t
% as the linear function on the domain.
eval([indVarNameSpace, '=xt;']);

% If we only have one variable appearing in allVarNames, the problem is a
% scalar problem.
scalarProblem = length(allVarNames) == 1;

if ( iscellstr(initInput) )
    % INITINPUT will be a cellstr in the case of coupled systems, where the
    % initial condition is specified over multiple lines.
    
    % Initialization:
    % The order in which the unknown variables appear in the initial condition.
    order = [];
    % The corresponding conditions:
    guesses = [];
    
    % Match LHS of = with variables in allVarNames
    for initCounter = 1:length(initInput)
        currStr = initInput{initCounter};
        equalSign = find(currStr == '=');
        currVar = strtrim(currStr(1:equalSign-1));
        match = find(ismember(allVarNames, currVar) == 1);
        order = [order ; match];
        currGuess = strtrim(currStr(equalSign+1:end));
        guesses = [guesses ; {currGuess}];
    end
    
    % If order is still empty, that means that initial guess were
    % passed on the form '2*x', rather than 'u = 2*x'. Allow that for
    % scalar problems, throw an error otherwise.
    if ( isempty(order) && scalarProblem )
        guess = eval(vectorize(initInput{1}));
        % If we have a scalar guess, convert to a chebfun
        if ( isnumeric(guess) )
            guess = 0*xt + guess;
        end
    elseif ( length(order) == length(guesses) )
        % We have a guess to match every dependent variable in the
        % problem.
        guess = cell(length(order),1);
        for guessCounter = 1:length(guesses)
            guessLoc = find(order == guessCounter);
            tempGuess = eval(vectorize(guesses{guessLoc}));
            if ( isnumeric(tempGuess) )
                tempGuess = 0*xt + tempGuess;
            end
            guess{guessCounter} = tempGuess;
        end
        
        % Convert the initial conditions to a CHEBMATRIX:
        guess = chebmatrix(guess);
    else
        % Something has gone wrong...
        error('CHEBFUN:CHEBGUI:constructInit:initMismatch', ...
            ['Error constructing initial guess.  Please make sure ' ...
            'guesses are of the form u = 2*x, v = sin(x), ...']);
    end
else
    % initInput is not a cell string, must have only received one line.
    guessInput = vectorize(initInput);
    equalSign = find(guessInput == '=');
    if ( isempty(equalSign) )
        equalSign = 0;
    end
    guessInput = guessInput(equalSign+1:end);
    guess =  chebfun(guessInput, [a b]);
end

end
