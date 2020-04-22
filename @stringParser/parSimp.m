function str = parSimp(str)
%PARSIMP   Remove unnecessary parentheses from string inputs.
%   STROUT = PARSIMP(STRIN) returns the string STROUT, obtained by doing some
%   basic simplifications of the string STRIN, attempting to remove unnecessary
%   parenthesis, zeros, and consecutive +/- pairs.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% We can't really expect to simplify a one character string...
if ( length(str) < 2 )
    return
end

% Simplify until we can't simplify more!
oldlength = 0;
while ( length(str) ~= oldlength )
    oldlength = length(str);
    str = parSimpMain(str);
end

end

function str = parSimpMain(str)
%PARSIMMAIN   Do the actual simplification

% Find all locations of ( and ):
leftParLoc = strfind(str, '(');
rightParLoc = strfind(str, ')');

% Find weighted vector of operators:
opVec = (str == '-') + (str == '+') + 2*((str == '*') + (str == '/')) + ...
    3*(str == '^');

% Operator locations
opLoc = find(opVec);

% Location of characters and numbers
charLoc = regexp(str, '[A-Za-z0-9@]');
ltgtLoc = regexp(str, '[\<\>]');

% Error if the number of ( and ) are not equal
if ( length(leftParLoc) ~= length(rightParLoc) )
    error('CHEBFUN:STRINGPARSER:parSimp:parenthSimplify', ...
        'Incorrect number of parenthesis.');
end

% Store number of parenthesis
numOfPars = length(leftParLoc);
% Pair them together
pairsLoc = zeros(numOfPars, 2);
pairsVec = zeros(1, length(str));
parCounter = 1;
leftStack = [];
for strIndex = 1:length(str)
    if ( str(strIndex) == '(' ) % Push to the stack
        leftStack = [leftStack, strIndex]; %#ok<AGROW>
    elseif ( str(strIndex) == ')' ) % Pop from the stack
        % Write information to pairs
        pairsLoc(parCounter, 1) = leftStack(end);
        pairsLoc(parCounter, 2) = strIndex;
        pairsVec([leftStack(end) strIndex]) = parCounter;
        % Update stack and counter
        leftStack(end) = [];
        parCounter = parCounter + 1;
    end
end

for k = 1:numel(ltgtLoc)
    idx = find((pairsLoc(:, 1) < ltgtLoc(k)) & (pairsLoc(:, 2) > ltgtLoc(k)));
    [ignored, idx2] = max(pairsLoc(idx, 1));
    charLoc = sort([charLoc pairsLoc(idx(idx2), :)]);
    pairsLoc(idx(idx2),:) = [];
    numOfPars = numOfPars - 1;
end

leftParsLoc = pairsLoc(:, 1);
rightParsLoc = pairsLoc(:, 2);

    function removeExtras(varargin)
        % Remove things other than parentheses, like consecutive +/- pairs
        % and hanging zeros. If there's an input argument, then we don't
        % attempt to remove zeros.
        
        % Remove leading + signs
        if ( strcmp(str(1), '+') )
            str(1) = []; 
            opVec(1) = [];
            leftParsLoc = leftParsLoc - 1;
            rightParsLoc = rightParsLoc - 1;
            opLoc = opLoc - 1;
            charLoc = charLoc - 1;
        end

        % Remove consecutive +/- pairs
        opVec1 = ( opVec == 1 );
        dOV = ~[abs(diff(opVec1)) 1] & opVec1; % Find adjacent opVec == 1's
        pm = find(dOV, 1); % Get the first
        if ( ~isempty(pm) ) % If +/- pair exists, remove the first occurence
            if ( str(pm) == str(pm+1) )
                str(pm) = '+';
            else
                str(pm) = '-';
            end
            str(pm + 1) = [];
            % Need to update indices and operators
            opVec(pm+1) = [];
            leftParsLoc(leftParsLoc > pm) = leftParsLoc(leftParsLoc > pm) - 1;
            rightParsLoc(rightParsLoc > pm) = ...
                rightParsLoc(rightParsLoc > pm) - 1;
            opLoc(opLoc > pm) = opLoc(opLoc > pm) - 1;
            charLoc(charLoc > pm) = charLoc(charLoc > pm) - 1;
            removeExtras(1) % recurse
        end

        % +/- pairs recurse locally
        if ( nargin > 0 )
            return
        end

    end

removeExtras()

% Loop over the pairs
for pIndex = 1:numOfPars
    
    removeExtras()
    
    pLeft = leftParsLoc(pIndex);
    pRight = rightParsLoc(pIndex);
    nextOpLeft = max(opLoc(opLoc < pLeft));
    nextOpRight = min(opLoc(opLoc > pRight));

    if ( nextOpLeft == 0 )
        nextOpLeft = [];
    end

    if ( nextOpRight == 0 )
        nextOpRight = [];
    end
    
    if ( any(charLoc == pLeft - 1) )
        % If the character in str next to the left of the left parenthesis is a
        % letter or a number (e.g. from sin or log2), we have a function from
        % which we cannot remove the () pair.
        continue
    end
    
    % Just mask the interior of these parentheses
    mask = zeros(1, length(str));
    mask(pLeft+1:pRight-1) = true;
    idx = find((leftParsLoc > pLeft) & (rightParsLoc < pRight));
    for k = 1:numel(idx)
        mask(leftParsLoc(idx(k)):rightParsLoc(idx(k))) = false;
    end
    interiorOpVec = opVec.*mask;
    interiorOpVec(interiorOpVec==0) = inf;
    minOpInside = min(interiorOpVec(pLeft:pRight));
    
    % If not, we perform a check to see whether we can remove them.
    if ( (isempty(nextOpLeft) || (minOpInside >= opVec(nextOpLeft))) && ...
            (isempty(nextOpRight) || minOpInside >= opVec(nextOpRight)) )
        
        if ( (minOpInside == 1) && ~isempty(nextOpLeft) && ...
                strcmp(str(pLeft - 1), '-') )
            % Deal with the case of minus (switch interior pluses and minuses
%            continue % uncomment if we don't want to simplify -(a+b) = -a-b
            interiorPM = interiorOpVec == 1;
            interiorPM(pLeft+1) = 0;
            m = (str == '-') & interiorPM;
            p = (str == '+') & interiorPM;
            % But we have to be careful that we don't change 1e-2 to 1e+2!
            mPos = find(m);
            pPos = find(p);
            mThrow = str(mPos-1) == 'e';
            pThrow = str(pPos-1) == 'e';
            m(mPos(mThrow)) = 0;
            p(pPos(pThrow)) = 0;
            str(m) = '+';
            str(p) = '-';
        end

        if ( (minOpInside == 2) && ~isempty(nextOpLeft) && ...
                strcmp(str(pLeft - 1), '/') )
            % Deal with division.  Unlike subtraction, we don't intend to
            % generally change a/(b/c) into a*c/b, so just move on.
            % (See GitHub issue #2357.)
            continue;
        end

        % Remove parenthesis pairs
        str(pLeft) = [];
        str(pRight-1) = [];

        % Update indices and operators
        opVec(pLeft) = [];
        opVec(pRight-1) = [];
        leftParsLoc(pIndex) = NaN;
        rightParsLoc(pIndex) = NaN;
        leftParsLoc = updateLoc(leftParsLoc, pLeft, pRight);
        rightParsLoc = updateLoc(rightParsLoc, pLeft, pRight);
        opLoc = updateLoc(opLoc, pLeft, pRight);
        charLoc = updateLoc(charLoc, pLeft, pRight);
    end
    
end

% Remove leading + signs
if ( strcmp(str(1), '+') )
    str(1) = []; 
end

end

function out = updateLoc(in, pLeft, pRight)
out = in - (in > pLeft) - (in > pRight);
end
