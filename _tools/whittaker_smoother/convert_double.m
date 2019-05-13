function [varargout]=convert_double(varargin)
% % Keywords: convert, double, data type
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % Syntax: [varargout]=convert_double(varargin);
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % Description
% % 
% % This program converts various data types to double precision numeric 
% % data type.  The number of input variables and output variables is
% % arbitrary.  This program is especially useful when processing several
% % variables of compressed data stored as single precision data type.  
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Input Variables
% % 
% % Can be any set of comma delimited numeric variables.
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % Output Variables 
% % 
% % Must have the same number of numeric output variables as input
% % variables.  The number of output variables is limited to the number of
% % input variables.  
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% 
% Example='';
%
% % Example 1 two variables
%
% SP=single(randn(1, 50000));
% a=single(1.698);
% [SP, a]=convert_double(SP, a);
%
% % Example 2 several variables
%
% SP=single(randn(100, 5000));
% a=single(1.698);
% b=single(rand(100, 100));
% c=single(rand(1000, 100));
% d=single(rand(10000, 100));
% [SP, a, b, c, d]=convert_double(SP, a, b, c, d);
%
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% % Program was written by Edward L. Zechmann
% %      date 17 August   2007
% %  modified 13 August   2008     updated comments
% %
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % Please feel free to modify this code.
% % 
% % See also: double, single, 
% % 


for e1=1:nargin;

    ttype=class(varargin{e1});
    
    if ~isequal(ttype, 'double') && logical( e1 >= nargout )
        varargout{e1}=double(varargin{e1});
    else
        varargout{e1}=varargin{e1};
    end
    
end
