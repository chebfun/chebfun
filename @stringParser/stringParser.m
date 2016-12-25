classdef stringParser
%STRINGPARSER   Parse strings for CHEBGUI.
%   This class is not intended to be called directly by the end user.
%
% See also CHEBGUI.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%   The @STRINGPARSER class implements a number of methods, used to convert a
%   string on a 'natural syntax' format in the fields of CHEBGUI to a format
%   that Chebfun is capable of working with. In v4, this functionality used to
%   live in the @chebgui folder, but to increase modularity, it has be spun off
%   to its own class.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function A = stringParser(varargin)
            % We never construct objects of this type, so the constructor is 
            % trivial.
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        % Convert a string to an anonymous function.
        varargout = str2anon(str, problemType, fieldType)
             
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = private, Static = true )
        
        % The STRINGPARSER lexer
        [out, varNames, pdeVarNames, eigVarNames, indVarNames] = ...
            lexer(str, type)
        
        % The LL(1) parser
        parseOut = parser(lexIn)
        
        % Convert an expression on PREFIX form to INFIX form
        [infixOut, notaVAR] = pref2inf(prefIn)
        
        % Split expressions separated by commas to individual syntax trees
        treeOut = splitTree(treeIn)
        
        % Split an EIG syntax tree
        [newTree, lambdaTree, lambdaSign] = splitTreeEIG(treeIn)
        
        % Split a PDE syntax tree
        [newTree, pdeSign] = splitTreePDE(treeIn)
        
        % Convert a syntax tree to prefix format
        prefixOut = tree2prefix(syntaxTree)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HIDDEN STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = public, Static = true, Hidden = true )
        
        % Get rid of unnecessary parenthesis in infix format strings
        strOut = parSimp(strIn)
        
    end
    
end
