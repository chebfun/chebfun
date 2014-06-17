function f = scribble(s, dom)
%SCRIBBLE   Write text with a complex-valued CHEBFUN.
%  F = JSCRIBBLE('STRING') returns a complex CHEBFUN representing the text in
%  STRING. F is paramterized by a real variable on the interval [-1, 1] and
%  maps the range of F is contained in the box [-1 1]x[0 1] in the complex
%  plane.
%
%  F = SCRIBBLE('STRING', DOM) paramterises F by a real variable on the
%  interval DOM.
%
%  The full United States QWERTY keyboard layout is supported (though
%  letters will be printed only in upper-case), as well as some special and
%  international characters. To contribute a character, please write it and
%  email it to help@chebfun.org.
%
%  Example:
%   f = scribble('The quick brown fox jumps over the lazy dog. 0123456789');
%   plot(f), axis equal

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

lengthS = length(s);
fData = {};
ns = 0; % Store the total number of strokes.
for j = 1:lengthS
    switch upper(s(j))
        % Alphabet
        case {'A'}, t = c([0 .4+1i .8 .6+.5i .2+.5i]);
        case {'B'}, t = c([0 1i .8+.9i .8+.6i .5i .8+.4i .8+.1i 0]);
        case {'C'}, t = c([.8+1i .8i .2i .8]);
        case {'D'}, t = c([0 .8+.1i .8+.9i 1i 0]);
        case {'E'}, t = c([.8+1i 1i .5i .5i+.7 .5i 0 .8]);
        case {'F'}, t = c([.8+1i 1i .5i .5i+.7 .5i 0]);
        case {'G'}, t = c([.8+1i .8i .2i .6 .6+.5i .4+.5i .8+.5i]);
        case {'H'}, t = c([0 1i .5i .5i+.8 .8+1i .8]);
        case {'I'}, t = c([0 .8 .4 .4+1i 1i .8+1i]);
        case {'J'}, t = c([0 .4 .4+1i 1i .8+1i]);
        case {'K'}, t = c([0 1i .5i .8+1i .5i .8]);
        case {'L'}, t = c([1i 0 .8]);
        case {'M'}, t = c([0 .1+1i .4 .7+1i .8]);
        case {'N'}, t = c([0 1i .8 .8+1i]);
        case {'O'}, t = c([0 1i .8+1i .8 0]);
        case {'Q'}, t = c([0 1i .8+1i .8 .6+.2i .9-.1i .8 0]);
        case {'P'}, t = c([0 1i .8+1i .8+.5i .5i]);
        case {'R'}, t = c([0 1i .8+1i .8+.6i .5i .8]);
        case {'S'}, t = c([.8+1i .9i .6i .8+.4i .8+.1i 0]);
        case {'T'}, t = c([.4 .4+1i 1i .8+1i]);
        case {'U'}, t = c([1i .1 .7 .8+1i]);
        case {'V'}, t = c([1i .4 .8+1i]);
        case {'W'}, t = c([1i .2 .4+1i .6 .8+1i]);
        case {'X'}, t = c([1i .8 .4+.5i .8+1i 0]);
        case {'Y'}, t = c([1i .4+.5i .8+1i .4+.5i .4]);
        case {'Z'}, t = c([1i .8+1i 0 .8]);
            % Numerals
        case {'0'}, t = c([0 1i .8+1i 0 .8 .8+1i]);
        case {'1'}, t = c([0 .8 .4 .4+1i .1+.8i]);
        case {'2'}, t = c([.8 0 .5i .8+.5i .8+1i 1i]);
        case {'3'}, t = c([1i .8+1i .8+.5i .1+.5i .8+.5i .8 0]);
        case {'4'}, t = c([1i .5i .8+.5i .8+1i .8]);
        case {'5'}, t = c([.8+1i 1i .5i .8+.5i .8 0]);
        case {'6'}, t = c([.8+1i 1i 0 .8 .8+.5i .5i]);
        case {'7'}, t = c([1i .8+1i .2]);
        case {'8'}, t = c([1i .8+1i .8 0 .5i 1i .5i .8+.5i]);
        case {'9'}, t = c([.8 .8+1i 1i .5i .8+.5i]);
            % Punctuation
        case {'.'}, t = c([0 .05 .05+.05i .05i 0]);
        case {','}, t = c([-0.1-.15i -.05-.15i .1+.05i .05+.05i -.1-.15i]);
        case {';'}, t = [c([-0.1-.15i -.05-.15i .1+.05i .05+.05i -.1-.15i]+.375+.3i) c([.05 .1 .1+.05i .05+.05i .05]+.375+.7i)];
        case {':'}, t = [c([0 .05 .05+.05i .05i 0]+.375+.3i) c([0 .05 .05+.05i .05i 0]+.375+.7i)];
        case {'?'}, t = [c([.6i 1i .8+1i .8+.5i .4+.5i .4+.2i]) c([.35 .45 .45+.05i .35+.05i .35])];
        case {'!'}, t = [c([.4+1i .4+.2i]) c([.35 .45 .45+.05i .35+.05i .35])];
        case {''''}, t = c([.3+.7i .4+1i]);
        case {'"'}, t = [c([.33+.8i .33+1i]) c([.66+.8i .66+1i])];
        case {'`'}, t = c([.5+.7i .4+1i]);
        case {'_'}, t = c([-.1 .9]);
        case {' '}, t = []; % Space
            % Parentheses
        case {')'}, t = c([.3+1i .4+.9i .5+.75i .5+.25i .4+.1i .3]);
        case {'('}, t = c([.5+1i .4+.9i .3+.75i .3+.25i .4+.1i .5]);
        case {'['}, t = c([.5 .3 .3+1i .5+1i]);
        case {']'}, t = c([.3 .5 .5+1i .3+1i]);
        case {'{'}, t = c([.5 .4+.1i .4+.4i .3+.5i .4+.6i .4+.9i .5+1i]);
        case {'}'}, t = c([.3 .4+.1i .4+.4i .5+.5i .4+.6i .4+.9i .3+1i]);
            % Symbols (math)
        case {'-'}, t = c([.1+.5i .7+.5i]);
        case {'+'}, t = c([.5i .8+.5i .4+.5i .4+.9i .4+.1i]);
        case {'*'}, t = [c([.4+.7i .4+.9i]) c([.32+.73i .48+.87i]) c([.48+.73i .32+.87i])];
        case {'/'}, t = c([.2 .6+1i]);
        case {'^'}, t = c([.2+.7i .4+1i .6+.7i]);
        case {'='}, t = [c([.1+.4i .7+.4i]) c([.1+.6i .7+.6i])];
        case {'<'}, t = c([.6+.3i .2+.5i .6+.7i]);
        case {'>'}, t = c([.2+.3i .6+.5i .2+.7i]);
            % Symbols (other)
        case {'\'}, t = c([.6 .2+1i]);
        case {'|'}, t = c([.4 .4+1i]);
        case {'%'}, t = [c([.1 .7+1i]) c(3*[0 .05 .05+.05i .05i 0]+.2+.85i) c(3*[0 .05 .05+.05i .05i 0]+.5+.05i)];
        case {'#'}, t = [c([.3+.1i .3+.9i]) c([.5+.1i .5+.9i]) c([.2+.3i .6+.3i]) c([.2+.7i .6+.7i])];
        case {'~'}, t = c([.2+.45i .3+.55i .5+.45i .6+.55i]);
        case {'@'}, t = c([.8 0 1i .8+1i .8+.25i .2+.25i .2+.75i .6+.75i .6+.25i]);
        case {'$'}, t = [c([.8+1i .9i .6i .8+.4i .8+.1i 0]) c([.4 .4+1i])];
        case {'£', '�'}, t = c([.8+1i .3+1i .2+.5i .5+.5i .5i .2+.5i .2 0.1 .8]);
        case {'&'}, t = c([.8 .7i 1i .7+1i .7+.7i .3i 0 .4 .8+.27i]);
            
        otherwise,
            t = [];
            warning('CHEBFUN:scribble:unknownChar', ...
                ['"', s(j), '" is not supported by scribble. (But please feel ', ...
                'free write it and email it to help@chebfun.org)']);
    end
    
    if ( ~isempty(t) )
        fData = [fData , t];
    end
    
end

% Construct the CHEBFUN representation of the text:
f = chebfun(fData, (0:ns));

% Rescale the ends to [-1 1] (or a given domain):
if ( nargin < 2 )
    dom = [-1, 1];
end
f = newDomain(f, dom);

    function cv = c(v)
        cv = cell(1, length(v)-1);
        for l = 2:length(v)
            % Shift:
            stroke = (v(l-1:l).' + (j-1))*2/lengthS - 1;           
            % Grab the x- and y-coordinates:
            sr = real(stroke);
            si = imag(stroke);
            % Construct the function handle:
            cv{l-1} = @(t)(t-ns)*(sr(2)-sr(1))+sr(1)+1i*((t-ns)*(si(2)-si(1))+si(1));
            % Count the total number of the strokes:
            ns = ns + 1;            
        end
    end

end
