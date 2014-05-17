classdef stringParser
%STRINGPARSER   Parse strings for CHEBGUI.
%   This class is not intended to be called directly by the end user.
%
%   See also CHEBGUI.
    
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    methods (Access = public)
        
        function A = stringParser(varargin)
            % We never construct objects of this type, so the constructor is
            % trivial
        end
        
    end
    
    methods( Access = public, Static = true )
        
        % Convert a string to an anonymous function
        varargout = str2anon(guifile, str, type)
        
        % The STRINGPARSER lexer
        [out, varNames, pdeVarNames, eigVarNames, indVarNames] = ...
            lexer(str, type)
    
        % The LL(1) parser
        parseOut = parser(lexIn)
        
        % Get rid of unecessary parenthesis in infix format strings
        strOut = parSimp(strIn)
        
        % Convert an expression on PREFIX form to INFIX form
        [infixOut, notaVAR] = pref2inf(prefIn)
        
        % Convert a syntax tree to prefix format
        prefixOut = tree2prefix(syntaxTree)
        
        % Split expressions separated by commas to individual syntax trees
        treeOut = splitTree_commas_equalSigns(treeIn)
        
        % Split a BVP syntax tree
        treeOut = splitTreeBVP(treeIn)
        
        % Split an EIG syntax tree
        [newTree, lambdaTree, lambdaSign] = splitTreeEIG(treeIn)
        
        % Split a PDE syntax tree
        [newTree, pdeSign] = splitTreePDE(treeIn)
    end
end