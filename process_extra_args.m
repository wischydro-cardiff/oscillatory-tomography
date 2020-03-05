function [varargout] = process_extra_args(orig_varargin,varargin)

%process_extra_args: Function used to take a varargin and either keep the
%default values for parameters or over-write them with user input from
%varargin.
%
% © 2005-2013 Michael Cardiff, All Rights Reserved.
%
%This function is a shortcut for writing a set of optional arguments for a
%function to their values. For example, if a function is defined as a =
%function_name(b,[c],[d]) where c and d are optional inputs, this can be
%dealt with as follows.
%
%Step 1) Create the default values for c and d within the program
%function_name
%Step 2) Within function_name, call process_extra_args as follows: [c d] =
%process_extra_args(varargin,c,d)
%
% process_extra_args will retain the default values for c and d unless
% there are non-blank elements in varargin, in which case these values will
% over-write the defaults.
%
% Syntax:
% [a b c d ...] = process_extra_args(orig_varargin,adef,bdef,cdef,ddef,...)
% where:
%   -a,b,c,d,... are the optional parameters for a function
%   -orig_varargin is the varargin array that has been passed to the
%   function
%   -adef,bdef,cdef,ddef,... are the default values for a,b,c,d,...
%   respectively.

nargs_toassign = nargin - 1;

nargs_supplied = numel(orig_varargin);

varargout = varargin;

for i = 1:1:nargs_toassign
    if i <= nargs_supplied
        if ~isempty(orig_varargin{i})
            varargout{i} = orig_varargin{i};
        end
    end
end