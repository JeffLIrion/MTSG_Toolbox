function display(GP)
% Display info about the GraphPart ojbect
%
% Input
%   GP      a GraphPart object
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% Input name
fprintf('\n\nGraphPart object ''%s''\n\n\n',inputname(1));


%% points & levels
[N, jmax] = size(GP.rs);
N = N-1;

if N > 0 && jmax > 0
    points = sprintf('|      # Points:  %d',N);
    spaces = floor(log10(N))-floor(log10(jmax));
    levels = sprintf('|      # Levels:  %s%d  <==>  jmax = %d',blanks(spaces),jmax,jmax-1);
end


%% class(rs) and class(ind)
rsind = sprintf('|      rs & ind:  %s',class(GP.rs));


%% method
method = sprintf('|        Method:  %s',GP.method);


%% tag
if ~isempty(GP.tag)
    tag = sprintf('|         "tag"?  Yes (%s)',class(GP.tag));
else
    tag = '|         "tag"?  No';
end



%% compinfo
if ~isempty(GP.compinfo)
    compinfo = '|    "compinfo"?  Yes';
else
    compinfo = '|    "compinfo"?  No';
end


%% rsf2c
if ~isempty(GP.rsf2c)
    rsf2c = '|       "rsf2c"?  Yes';
else
    rsf2c = '|       "rsf2c"?  No';
end


%% tagf2c
if ~isempty(GP.tagf2c)
    tagf2c = '|      "tagf2c"?  Yes';
else
    tagf2c = '|      "tagf2c"?  No';
end


%% compinfof2c
if ~isempty(GP.compinfof2c)
    compinfof2c = '| "compinfof2c"?  Yes';
else
    compinfof2c = '| "compinfof2c"?  No';
end


%% print the info
if exist('points','var')
    tablewidth = max([length(points), length(levels), length(rsind), length(method), length(tag), length(compinfo), length(rsf2c), length(tagf2c), length(compinfof2c)]) + 2;
    
    fprintf(' %s \n',repmat('-',1,tablewidth-2));
    fprintf('|%s|\n',blanks(tablewidth-2));
    fprintf('%s%s|\n',points,blanks(tablewidth-length(points)-1));
    fprintf('%s%s|\n',levels,blanks(tablewidth-length(levels)-1));
    fprintf('%s%s|\n',rsind,blanks(tablewidth-length(rsind)-1));
    fprintf('|%s|\n',blanks(tablewidth-2));
    fprintf('%s%s|\n',method,blanks(tablewidth-length(method)-1));
    fprintf('|%s|\n',blanks(tablewidth-2));
    fprintf('%s%s|\n',tag,blanks(tablewidth-length(tag)-1));
    fprintf('%s%s|\n',compinfo,blanks(tablewidth-length(compinfo)-1));
    fprintf('|%s|\n',blanks(tablewidth-2));
    fprintf('%s%s|\n',rsf2c,blanks(tablewidth-length(rsf2c)-1));
    fprintf('%s%s|\n',tagf2c,blanks(tablewidth-length(tagf2c)-1));
    fprintf('%s%s|\n',compinfof2c,blanks(tablewidth-length(compinfof2c)-1));
    fprintf('|%s|\n',blanks(tablewidth-2));
    fprintf(' %s ',repmat('-',1,tablewidth-2));
    
end
    
fprintf('\n\n\n');


end