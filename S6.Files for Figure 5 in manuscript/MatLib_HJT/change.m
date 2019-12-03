function new = change(old,relation,flag,value)

% CHANGE     Change values in a matrix
%--------------------------------------------------------------------
% CHANGE   1.3   92/03/25
%
% new = change(old,relation,flag,value)
%
% DESCRIPTION:
%    Changes the 'flag' values in the matrix "old" to the new "value"
%    according to the "relation".
% 
% INPUT:
%    old      = matrix containing values related to "flag"
%               are to be converted to "value"
%    flag     = values related to "flag" then replaced by "value"
%    value    = replacement value
%    relation = string relation e.g. '<', '>', '=='
%
% OUTPUT:
%    new      = matrix "old" with all flagged values changed
%               to "value" (can be returned to same matrix "old")
%
% EXAMPLE:  A = change(A,'==',NaN,0  )
%           B = change(A,'<', -99,Nan) 
%
% CALLER:   general purpose
% CALLEE:   none
%
% AUTHOR:   Phil Morgan 3-02-92

% @(#)change.m   1.3   92/03/25
% 
% Re-created after 2-2-92 hard disk crash.
% Based on flagnan.m - Steve Rintoul  Dec 90
%          alter.m   - Phil  Morgan    Feb 91
%          convert.m - Peter Mcintosh Aug 91
%
% John Wilkin added handling of ~= for flag=NaN and removed restriction on
% relation operators 
% 
%--------------------------------------------------------------------

% CHECK INPUT ARGUMENTS CALL
if nargin ~= 4
  error('CHANGE.M: Must have 4 input arguments')
end

if strcmp(relation,'==') | strcmp(relation,'>') | strcmp(relation,'<') | ...
   strcmp(relation,'~=') | strcmp(relation,'>=') | strcmp(relation,'<=')
    % valid relation
  else  
    error(['CHANGE.M: Relation {' relation '} not valid'])
end 

if isnan(flag)
   switch relation
      case '=='
        replace = find(isnan(old));
      case '~='
        replace = find(~isnan(old));
      otherwise
        error('CHANGE.M: Relation should be == or ~= to compare to NaN')
   end
else
   evalstr = ['replace = find(old',relation,'flag);'];
   % disp(evalstr)
   eval(evalstr)
end
new = old;
if isempty(value)
   new(replace) = [];
else
   new(replace) = repmat(value,size(replace));
end
