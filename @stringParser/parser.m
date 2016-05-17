function parseOut = parser(lexIn)
%STRCONVPARSER    LL(1) parser for mathematical expressions
%  PARSEOUT = STRCONVPARSER(LEXIN) returns a syntax tree of expressions so that
%  it can be converted to a format Chebfun is able to work with. The input,
%  LEXIN, is the output of the method STRCONVLEXER(), and is a cell array of
%  strings, containaining the tokens of strings.
%
%  This method implements a LL(1) parser, a technique from compiler theory. For
%  more details, see e.g.,
%   [1] Aho, Sethi, Ullman, Compilers: Principles, Techniques, and Tools,
%       Addison-Wesley, 1986.
%
% See also: STRINGPARSER, STRINGPARSER/STR2ANON, STRINGPARSER/LEXER.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developers note: Usually, a parser relies on having access to pointers, which
% is not possible in MATLAB. We get around this issue using global variables.
%
% TODO: Can we now get around this? See #515.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize all global variables
global NEXT
global COUNTER
global LEX
global STACK

% Enter the main routine
parseMain(lexIn);

% Return the stored stack
parseOut = STACK;

% Clear all global variables
NEXT = [];
COUNTER = [];
LEX = [];
STACK = [];

end

function parseMain(lexIn)
%PARSEMAIN  The main parsing routine, starts the recursive parsing.

global NEXT
global STACK
global LEX
global COUNTER

COUNTER = 1;
LEX = lexIn;
NEXT = char(LEX(COUNTER, 2));
STACK = [];

% Our expression can only start with certain labels, make sure we are
% starting with one of them.
validTypes = {'NUM', 'VAR', 'INDVAR', 'PDEVAR', 'LAMBDA',...
            'FUNC1', 'FUNC2', 'FUNC3', 'UN-', 'UN+', 'LPAR'};

if ( any(strcmp(NEXT, validTypes)) )
    
    % Enter the recursion and check for successful termination.
    parseExpA();
    
    % We've put a $ at the end of the lexer output. Hence, we have only
    % successfully parsed our string if the only remaining token left is a
    % $ sign.
    success = match('$');
    
    if ( ~success )
        reportError('Parse:end', ...
            'Input expression ended in unexpected manner.');
    end
    
else
    
    reportError('Parse:start', 'Input field started with unaccepted symbol.');
    
end

end

% The following is the allowed grammar of mathematical expressions in the
% chebgui (¬ denotes an empty symbol):
%
%   ExpA    -> ExpA   COMMA  Exp0
%   ExpA    -> Exp0
%   Exp0    -> Exp0   OP=    Exp05
%   Exp0    -> Exp05
%   Exp05   -> Exp05  OPREL  Exp1       where OPREL can be OP>, OP>=, OP<, OP<=
%   Exp05   -> Exp1
%   Exp1    -> Exp1   OP+    Exp2
%   Exp1    -> Exp1   OP-    Exp2
%   Exp1    -> Exp2
%   Exp2    -> Exp2   OP*    Exp3
%   Exp2    -> Exp2   OP/    Exp3
%   Exp2    -> Exp3
%   Exp3    -> Exp3   OP^    Exp4
%   Exp3    -> Exp4
%   Exp4    ->        UN+    Exp4
%   Exp4    ->        UN-    Exp4
%   Exp4    -> Exp5
%   Exp5    -> NUM
%   Exp5    -> INDVAR
%   Exp5    -> PDEVAR
%   Exp5    -> LAMBDA
%   Exp5    -> VAR
%   Exp5    -> VAR ExpDer
%   Exp5    -> VAR LPAR ExpList RPAR
%   Exp5    -> VAR ExpDer LPAR ExpList RPAR
%   Exp5    -> FUNC1 LPAR Exp1 RPAR
%   Exp5    -> FUNC2 LPAR Exp1 COMMA Exp1 RPAR]
%   Exp5    -> FUNC3 LPAR Exp1 COMMA Exp1 COMMA Exp1 RPAR
%   Exp5    -> LPAR Exp1 RPAR
%   ExpDer  -> DER
%   ExpDer  -> ¬
%   ExpList -> Exp1
%   ExpList -> Exp1 COMMA STR

function parseExpA()

parseExp0();
parseExpApr();

end

function parseExp0()

parseExp05();
parseExp0pr();

end

function parseExp05()

parseExp1();
parseExp05pr();

end

function parseExp1()

parseExp2();
parseExp1pr();

end

function parseExp2()

parseExp3();
parseExp2pr();

end

function parseExp3()

parseExp4();
parseExp3pr();

end

function parseExp4()

parseExp5();

end

function parseExp5()

global NEXT
global COUNTER
global LEX

% We begin by checking whether we have hit a terminal case. In that case,
% we push that into the stack. We need to treat variables with _ in the
% names separately, as we only allow certain operators around the time
% derivative.
if ( strcmp(NEXT, 'VAR') )
    % Create a new tree for later use, corresponding to the dependent
    % variable
    tempLeaf = struct('center', {{char(LEX(COUNTER)), char(NEXT)}}, ...
        'pdeflag', 0);
    
    % Push the variable onto the stack so that we can take the correct
    % actions if we have derivatives or expressions on the form u(...)
    % involved. If not, we'll simply pop it later.
    push(tempLeaf);

    % Begin by advancing as usual.
    advance();
    
    % If there follows a derivative symbol, ', we parse that
    parseExpDer();
    
    % If we then have a u(3) or u(3,left) situation (or u'(3) etc.), we
    % parse that as required
    parseExpList();
    
elseif ( any(strcmp(NEXT, {'NUM', 'INDVAR', 'LAMBDA', 'STR'})) )
    newLeaf = struct('center', {{char(LEX(COUNTER)), char(NEXT)}}, ...
        'pdeflag', 0);
    push(newLeaf);
    advance();
    
elseif ( strcmp(NEXT, 'PDEVAR') )
    newLeaf = struct('center', {{char(LEX(COUNTER)), char(NEXT)}}, ...
        'pdeflag', 1);
    push(newLeaf);
    advance();
    
elseif ( strcmp(NEXT, 'FUNC1') ) % Functions which take one argument
    parseFunction1();
    
elseif ( strcmp(NEXT, 'FUNC2') ) % Functions which take two arguments
    parseFunction2();
    
elseif ( strcmp(NEXT, 'FUNC3') ) % Functions which take three arguments
    parseFunction3();

elseif ( strcmp(NEXT, 'LPAR') )
    advance();
    parseExp05();

    % Check if NEXT symbol is ')' as it should be. If not, there is a
    % parenthesis imbalance in the expression and we return an error.
    m = match('RPAR');  
    if ( ~m )
        reportError('Parse:parenths', 'Parenthesis imbalance in input fields.')
    end
    
elseif  ( strcmp(NEXT, 'UN-') || strcmp(NEXT, 'UN+') || ...
          strcmp(NEXT, 'OP-') || strcmp(NEXT,'OP+') )
    % If + or - reaches this far, we have an unary operator.
    % ['UN', char(NEXT(3))] determines whether we have UN+ or UN-.
    newCenterNode = {{char(LEX(COUNTER)), ['UN', char(NEXT(3))]}};
    advance();
    parseExp4();
    
    rightArg = pop();
    
    pdeflag = rightArg.pdeflag;

    newTree = struct('center', newCenterNode, 'right', rightArg, ...
        'pdeflag',pdeflag);
    push(newTree);
else
    
    reportError('Parse:terminal', ...
        ['Unrecognized character in input field:', NEXT]);
end

end

function parseFunction1()

global COUNTER
global LEX

functionName = char(LEX(COUNTER));
advance();

if ( match('LPAR') )
    parseExp1();
    
    if ( match('COMMA') )
        reportError('Parse:func1', ...
            ['Method ''', functionName, ''' only takes one input argument.']);
    elseif ( ~match('RPAR') )
        reportError('Parse:parenths', 'Parenthesis imbalance in input fields.')
    end
    
    rightArg =  pop();
    if ( rightArg.pdeflag )
        error('CHEBFUN:STRINGPARSER:parser:PDE', ...
            'Cannot use time derivative as function arguments.')
    end
    % Can assume no pde if we reach here
    newTree = struct('center', {{functionName, 'FUNC1'}}, ...
        'right', rightArg, 'pdeflag', 0);
    push(newTree);
else
    reportError('Parse:parenths', ...
        'Need parenthesis when using functions in input fields.')
end

end

function parseFunction2()

global NEXT
global COUNTER
global LEX

% Function which allow one or two arguments
oneArgAllowed = {'diff', 'cumsum', 'airy', 'mean'};
functionName = char(LEX(COUNTER));
advance();

% Here we need ( as the next symbol
if ( strcmp(NEXT, 'LPAR') )
    advance();
    parseExp1();
    
    firstArg =  pop();
    if ( firstArg.pdeflag )
        error('CHEBFUN:STRINGPARSER:parser:PDE', ...
            'Cannot use time derivative as function arguments.')
    end
    
    % Check whether we have a comma, if so, continue as normal
    if ( match('COMMA') )
        parseExp1();
        m = match('RPAR');
        if ( ~m )
            reportError('Parse:parenths', ...
                'Parenthesis imbalance in input fields.')
        end
        
        secondArg =  pop();
        if ( secondArg.pdeflag )
            error('CHEBFUN:STRINGPARSER:parser:PDE', ...
                'Cannot use time derivative as function arguments.')
        end

        % Can assume no pde if we reach here
        newTree = struct('left', firstArg, ...
            'center', {{functionName, 'FUNC2'}}, ...
            'right', secondArg,'pdeflag',0);
        push(newTree);
        
    elseif ( match('RPAR') )
        if ( any(strcmp(functionName, oneArgAllowed)) )
            % Have hit a function which allows one or two args.

            % If we only had one argument, we convert the function to type
            % FUNC1.  Can assume no PDE if we reach here
            newTree = struct('center', {{functionName, 'FUNC1'}}, ...
                'right', firstArg, 'pdeflag', 0);
            push(newTree);
        else
            % We tried to call a method which requires two args with only one.
            reportError('Parse:func2', ['Method ''', functionName, ...
                ''' requires two input arguments.']);
        end
        
    else
        reportError('Parse:parenths', 'Parenthesis imbalance in input fields.')
    end
    
else
    reportError('Parse:parenths', ...
        'Need parenthesis when using functions in input fields.')
end

end

function parseFunction3()

global NEXT
global COUNTER
global LEX

% Function which allow one or two arguments
oneArgAllowed = {'sum', 'integral'};
twoArgAllowed = {'feval', 'fred', 'volt'};
functionName = char(LEX(COUNTER));
advance();

% Here we need ( as the next symbol
if ( strcmp(NEXT, 'LPAR') )
    advance();
    parseExp1();
    firstArg = pop();
    if ( firstArg.pdeflag )
        error('CHEBFUN:STRINGPARSER:parser:PDE', ...
            'Cannot use time derivative as function arguments.')
    end
    
    % Check whether we have a comma, if so, continue as normal
    if ( match('COMMA') )
        parseExp1();
        secondArg = pop();
        if ( secondArg.pdeflag )
            error('CHEBFUN:STRINGPARSER:parser:PDE', ...
                'Cannot use time derivative as function arguments.')
        end

        if ( match('COMMA') ) % and again
            parseExp1();
            thirdArg = pop();

            % Define the new branch
            newTree = struct('left', firstArg, ...
                'center', {{functionName, 'FUNC3'}}, ...
                'right', secondArg, 'arg', thirdArg, 'pdeflag', 0);

            if ( ~match('RPAR') ) % Check the final parenthesis
                reportError('Parse:parenths', ...
                    'Parenthesis imbalance in input fields.')
            end
        elseif ( match('RPAR') && any(strcmp(functionName, twoArgAllowed)) )
            % This was actually a two argument function in disguise
            newTree = struct('left', firstArg, ...
                'center', {{functionName, 'FUNC2'}},...
                'right', secondArg, 'pdeflag', 0);
        else
            reportError('Parse:parenths', ...
                'Parenthesis imbalance in input fields.')
        end
        push(newTree);
        
    elseif ( match('RPAR') )
        if ( any(strcmp(functionName,oneArgAllowed)) )
            % Have hit a function which allows one or two args.

            % If we only had one argument, convert the function to type FUNC1.
            % Can assume no PDE if we reach here.
            newTree = struct('center', {{functionName, 'FUNC1'}}, ...
                'right', firstArg, 'pdeflag', 0);
            push(newTree);
        else
            % We tried to call a method which requires two args with only one.
            reportError('Parse:func2', ['Method ''', functionName, ...
                ''' requires two input arguments.']);
        end
        
    else
        reportError('Parse:parenths', 'Parenthesis imbalance in input fields.')
    end
    
else
    reportError('Parse:parenths', 'Need parenthesis when using functions in input fields.')
end

end

function parseExp1pr()

global NEXT

if ( strcmp(NEXT, 'OP+') )

    advance();
    leftArg  = pop();
    parseExp2();
    rightArg = pop();
    
    pdeflag = leftArg.pdeflag || rightArg.pdeflag;
    
    newTree = struct('left', leftArg, 'center', {{'+', 'OP+'}}, ...
        'right', rightArg, 'pdeflag', pdeflag);
    push(newTree);
    parseExp1pr();
    
elseif ( strcmp(NEXT, 'OP-') )
    advance();
    leftArg  = pop();
    parseExp2();
    rightArg = pop();

    pdeflag = leftArg.pdeflag || rightArg.pdeflag;
    newTree = struct('left', leftArg, 'center', {{'-', 'OP-'}}, ...
        'right', rightArg, 'pdeflag', pdeflag);
    push(newTree);
    parseExp1pr();
elseif ( strcmp(NEXT, 'RPAR') || strcmp(NEXT, '$') || strcmp(NEXT, 'OP=') )
	% Do nothing
else % If we don't have ) or the end symbol now something has gone wrong.
%     reportError('Parse:end','Syntax error in input fields.')
end

end


function parseExp2pr()

global NEXT

if ( strcmp(NEXT,'OP*') )
    leftArg  = pop();   % Pop from the stack the left argument
    advance();          % Advance in the input
    parseExp3();
    rightArg = pop();  % Pop from the stack the right argument

    % Check whether we have _ variables
    if ( leftArg.pdeflag || rightArg.pdeflag )
        error('CHEBFUN:STRINGPARSER:parser:PDE','Cannot multiply time derivative')
    end

    % Can assume no PDE if we reach here
    newTree = struct('left', leftArg, 'center', {{'.*', 'OP*'}}, ...
        'right', rightArg, 'pdeflag', 0);
    push(newTree);
    parseExp2pr();
    
elseif ( strcmp(NEXT,'OP/') )
    leftArg  = pop();   % Pop from the stack the left argument
    advance();          % Advance in the input
    parseExp3();
    rightArg = pop();  % Pop from the stack the right argument

    % Check whether we have _ variables
    if ( leftArg.pdeflag || rightArg.pdeflag )
        error('CHEBFUN:STRINGPARSER:parser:PDE','Cannot divide with time derivatives')
    end

    % Can assume no PDE if we reach here
    newTree = struct('left', leftArg, 'center', {{'./', 'OP/'}}, ...
        'right', rightArg, 'pdeflag', 0);
    push(newTree);
    parseExp2pr();
else
    % Do nothing
end

end

function parseExp3pr()

global NEXT

if ( strcmp(NEXT,'OP^') )
    leftArg  = pop();
    advance();    
    parseExp4();
    rightArg = pop();

    % Check whether we have _ variables
    if ( leftArg.pdeflag || rightArg.pdeflag )
        error('CHEBFUN:STRINGPARSER:parser:PDE','Cannot take powers with time derivative')
    end

    % Can assume no pde if we reach here
    newTree = struct('left', leftArg, 'center', {{'.^', 'OP^'}}, ...
        'right', rightArg, 'pdeflag', 0);
    push(newTree);
    parseExp3pr();
    
elseif ( ~isempty(strfind(NEXT, 'DER')) )
    leftArg  = pop();

    % Check whether we have _ variables
    if ( leftArg.pdeflag )
        error('CHEBFUN:STRINGPARSER:parser:PDE','Cannot differentiate time derivative')
    end

    newTree = struct('center', {{'D', NEXT}}, 'right', leftArg, 'pdeflag', 0);
    push(newTree);
    advance();
    parseExp3pr();
else
    % Do nothing
end

end

function parseExpDer()
% parseExpDer deals with ' denoting derivatives on the dependent variables
% in the problem.

global NEXT

if ( ~isempty(strfind(NEXT, 'DER')) )
    leftArg  = pop();

    % Check whether we have _ variables
    if ( leftArg.pdeflag )
        error('CHEBFUN:STRINGPARSER:parser:PDE','Cannot differentiate time derivative')
    end

    newTree = struct('center', {{'D', NEXT}}, 'right', leftArg, 'pdeflag', 0);
    push(newTree);
    advance();
else
    % Do nothing
end

end

function parseExpList()
% parseExpList deals with potential arguments to the dependent variables in
% the problem, e.g. u(3,left).

global NEXT
global COUNTER
global LEX

if ( match('LPAR') ) % We are in a u(3) or u(3,left) situation
    % The first argument expression can only be of type EXP05 or higher
    % (no = or comma-dividers are allowed)
    parseExp1();
    % Check whether we are in in u(3) or u(3,left) situation. If we
    % have a match with ), we must have a feval with two arguments. If
    % we have a match with a comma, we must have a feval with three
    % arguments, but if there is no match, we have parenthesis
    % imbalance.
    if ( match('RPAR') )
        % Here we only had one argument in u(...), so we create a feval of type
        % FUNC2 (since the feval method has two input arguments)
        secondArg = pop();
        firstArg = pop();

        % Can assume no PDE if we reach here
        newTree = struct('left', firstArg, 'center', {{'feval', 'FUNC2'}}, ...
            'right', secondArg,'pdeflag', 0);
        push(newTree);
        
    elseif ( match('COMMA') )
        % Here we only had two arguments in u(...), so we create a
        % feval of type FUNC3 (since the feval method has three input
        % arguments)
        
        % Check whether we have any of the allowed option type as the
        % second argument. If not, throw an error.
        if ( strcmp(NEXT,'STR') )
            % We got a match! But we also need to check for parenthesis
            % balance. Store the option as a temporary leaf
            optLeaf = struct('center', {{char(LEX(COUNTER)), char(NEXT)}}, ...
                'pdeflag', 0);
            % Advance and check for parenth. balance
            advance();
            if ( match('RPAR') )
                % Create a feval of type FUNC3 (since the feval method
                % has three input arguments)
                secondArg = pop();
                firstArg = pop();
                thirdArg = optLeaf;

                % Can assume no pde if we reach here
                newTree = struct('left', firstArg, ...
                    'center', {{'feval', 'FUNC3'}}, 'right', secondArg, ...
                    'arg',thirdArg, 'pdeflag', 0);
                push(newTree);
            else
                reportError('Parse:parenths', ...
                    'Parenthesis imbalance in input fields.')
            end
        else
            reportError('Parse:secondArg', ...
                'Invalid second argument to u(0,...) type of expression.')
        end
        
    else % There must be a parenthesis imbalance.
        reportError('Parse:parenths', 'Parenthesis imbalance in input fields.')
    end
else % We reach here for the default behaviour, e.g. u' = 1, no parentheses involved
    % Do nothing, since we've already pushed the variable onto the
    % stack.
end

end

function parseExp0pr()

global NEXT

if ( strcmp(NEXT, 'OP=') )
    
    leftArg  = pop();
    advance();
    parseExp05();

    rightArg = pop();

    pdeflag = leftArg.pdeflag || rightArg.pdeflag;

    newTree = struct('left', leftArg, 'center', {{'=', 'OP='}}, ...
        'right', rightArg, 'pdeflag', pdeflag);
    push(newTree);   
else
	% Do nothing
end

end

function parseExpApr()

global NEXT

if ( strcmp(NEXT, 'COMMA') )
    
    advance();
    leftArg  = pop();
    parseExpA();
    rightArg = pop();
    
    pdeflag = leftArg.pdeflag || rightArg.pdeflag;
    
    newTree = struct('left', leftArg, 'center', {{',', 'COMMA'}}, ...
        'right', rightArg, 'pdeflag', pdeflag);
    push(newTree);
else
	% Do nothing
end

end

function parseExp05pr()

global NEXT

if ( any(strcmp(NEXT, {'OP>', 'OP>=', 'OP<', 'OP<='})) )
    tempOpType = NEXT;
    tempOpLabel = tempOpType(3:end);
    
    advance();
    leftArg  = pop();
    parseExp1();
    rightArg = pop();
    
    pdeflag = leftArg.pdeflag || rightArg.pdeflag;
    
    newTree = struct('left', leftArg, 'center', {{tempOpLabel, tempOpType}}, ...
        'right', rightArg, 'pdeflag', pdeflag);
    push(newTree);
else
	% Do nothing
end

end

function advance()
%ADVANCE    Move to the next token in the output from the lexer.
global NEXT
global COUNTER
global LEX

COUNTER = COUNTER + 1;
NEXT = char(LEX(COUNTER,2));

end



function m = match(label)

global NEXT
m = strcmp(label, NEXT);

% If we found a match and are not at the end of output, we want to advance
% to the NEXT symbol
if ( m && ~strcmp(label, '$') )
    advance();
end

end


function push(new)
%PUSH   Push a tree to the stack of syntax trees
global STACK

if ( ~stackRows() )
    STACK = new;
else
    % Ensure number of fields matches
    [STACK, new] = mergeFields(STACK, new);
    % Update the Stack
    STACK = [STACK ; new];
end

end


function p = pop()
%POP    Pop a tree from the stack.
global STACK

% Throw a sensible error if we have an empty stack
if ( isempty(STACK) )
    reportError('Parse:end', 'Syntax error in input expression (Empty stack)');
end

p = STACK(end,:);
STACK(end,:) = [];

end

function m = stackRows()

global STACK
[m, n] = size(STACK); %#ok<NASGU>

end

% TODO:  Can we get rid of this function entirely?
function reportError(id, msg)

error(['CHEBFUN:', id], msg);
% ME = MException(id,msg);
% throw(ME);

end


function [a, b] = mergeFields(a, b)

aFields = fieldnames(a); 
bFields = fieldnames(b);

for k = 1:numel(aFields)
    if ( ~isfield(b, aFields{k}) )
        b.(aFields{k}) = [];
    end
end

for k = 1:numel(bFields)
    if ( ~isfield(a, bFields{k}) )
        a.(bFields{k}) = [];
    end
end

end
