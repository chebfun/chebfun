function varargout = why(f, r)
%WHY   Provides a succinct answer to almost any Chebfun-related question in the
%      many languages of the friends of Chebfun.
%   WHY, by itself, provides a random answer.
%   WHY(N) provides the N-th answer.
%
%   For fun, try also
%       plot(why(chebfun)), axis equal

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = 12;  % Number of Chebfun languages
if ( nargin == 1 )
    r = ceil(n*rand);
end

if ( verLessThan('matlab', '7.9') )
    simpleText = true;
else
    simpleText = false;
end

switch r
    case 1,
        % English (or American-English!):
        % Nick Trefethen, Toby Driscoll, Nick Hale,
        % Sheehan Olver, Mark Richardson, Zachary Battles
        % Alex Townsend, James Lottes, Matthew Green

        s = 'Because Nick said so!';

    case 2,
        % Spanish:
        % Ricardo Pachon

        s = '!Porque Nick lo dice!';

    case 3,
        % Portuguese:
        % Rodrigo Platte

        s = 'Porque o Nick disse!';

    case 4,
        % Icelandic:
        % Asgeir Birkisson

        if ( ~simpleText )
            s = 'Því Nick mælti svo!';
        else
            s = 'Thvi Nick maelti svo!';
        end

    case 5,
        % Afrikaans:
        % Andre Weideman

        if ( ~simpleText )
            s = 'Omdat Nick so sê!';
        else
            s = 'Omdat Nick so se!';
        end

    case 6,
        % Dutch:
        % Joris Van Deun

        s = 'Omdat Nick het zegt!';

    case 7,
        % Swiss German:
        % Pedro Gonnet

        s = 'Wills dr Nick gseit hett!';

    case 8,
        % French:
        % Cecile Piret

        s = 'Parce-que Nick le dit !';

    case 9,
        % Urdu:
        % Mohsin Javed

        s = 'Kiyon ke Nick ne kaha tha!';

    case 10,
        % Luxembourgish:
        % Georges Klein

        s = 'Well den Nick et gesot huet!';
        
    case 11,
        % German:
        % Stefan Guettel
        
        s = 'Weil es Nick so gesagt hat!';
        
    case 12,
        % Chinese:
        % Kuan Xu

        if ( ~simpleText )
            s = '因为尼克是这么说的！';
        else
            s = 'Yin Wei Ni Ke Shi Zhe Me Shuo De!';
        end

    otherwise,
        s = 'Good question!';

end

if ( nargout < 1 )
    disp(s)
else
    varargout{1} = scribble(s);
end
