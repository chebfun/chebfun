function varargout = why(f, r)  
%WHY   Provides a succinct answer to almost any Chebfun related question in the 
%      many languages of the friends of Chebfun.
%   WHY, by itself, provides a random answer.
%   WHY(N) provides the N-th answer.
%   
%   For fun, try also
%       plot(why(chebfun)), axis equal

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

N = 11;  % Number of Chebfun languages
if ( nargin == 1 )
    r = ceil(N*rand);
end

if ( verLessThan('matlab', '7.9') )
    simpleText = true;
else
    simpleText = false;
end
    
switch r
    case 1,         s = 'Because Nick Trefethen said so!';
        % English (or American-English!):
        % Nick Trefethen, Toby Driscoll, Nick Hale,  
        % Sheehan Olver, Mark Richardson, Zachary Battles
        % Alex Townsend, James Lottes, Matthew Green
        
    case 2,         s = 'Porque Nick Trefethen lo dice!';
        % Spanish:
        % Ricardo Pachon
        
    case 3,         s = 'Porque o Nick Trefethen disse!';
        % Portuguese:
        % Rodrigo Platte
                
    case 4,         
        if ( ~simpleText )
                    s = 'Því Nick Trefethen mælti svo!';
        else
                    s = 'Thvi Nick Trefethen maelti svo!';
        end
        % Icelandic:
        % Asgeir Birkisson
        
    case 5,
        if ( ~simpleText )
                    s = 'Omdat Nick Trefethen so sê!';
        else
                    s = 'Omdat Nick Trefethen so se!';
        end
        % Afrikaans:
        % Andre Weideman
        
    case 6,         s = 'Omdat Nick Trefethen het zegt!';
        % Dutch:
        % Joris Van Deun
        
    case 7,         s = 'Wills dr Nick Trefethen gseit hett!';
        % Swiss German:
        % Pedro Gonnet
    
    case 8,         s = 'Parce que Nick Trefethen le dit!';
        % French:
        % Cecile Piret
        
    case 9,         s = 'Kiyon ke Nick Trefethen ne kaha tha!'; 
        % Urdu:
        % Mohsin Javed
        
    case 10,        s = 'Well den Nick Trefethen et gesot huet!';
        % Luxembourgish:
        % Georges Klein
        
    case 11,        s = 'Weil es Nick Trefethen so gesagt hat!';
        % German:
        % Stefan Guettel

    otherwise,      s = 'Good question!';
        
end

if ( nargout < 1 )
    disp(s)
else
    varargout{1} = scribble(s);
end
