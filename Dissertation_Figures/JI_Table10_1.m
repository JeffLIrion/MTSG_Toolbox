% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



if ~exist('ScienceNews.mat','file')
    fprintf('\n\n\nIn order to generate these results, install Python,\nrun WavePath.m, then run this script again.\n\n\n');
    return
end

load('ScienceNews.mat');
[rows,cols] = size(cmatrix);

% tabulate the table
TABLE = tabulate(double(classes));

% edit the class names
class_names = [class_names, repmat(' ',8,22)];
for row = 1:8
    class_names(row,:) = strrep(class_names(row,:),'Anth                     ','Anthropology & Archeology');
    class_names(row,:) = strrep(class_names(row,:),'Space                     ','Astronomy & Space Sciences');
    class_names(row,:) = strrep(class_names(row,:),'Behav   ','Behavior');
    class_names(row,:) = strrep(class_names(row,:),'Env                           ','Earth & Environmental Sciences');
    class_names(row,:) = strrep(class_names(row,:),'Life         ','Life Sciences');
    class_names(row,:) = strrep(class_names(row,:),'MathCS                 ','Mathematics & Computers');
    class_names(row,:) = strrep(class_names(row,:),'Med             ','Medical Sciences');
    class_names(row,:) = strrep(class_names(row,:),'PhysTech                     ','Physical Science & Technology');
end


% print the table
fprintf('\n\n         Class                    # Documents    %% of total\n');
fprintf('%s\n',repmat('-',1,59));
for row = 1:length(TABLE)
    fprintf('%s       %4.0f           %5.2f%%\n',class_names(row,:),TABLE(row,2),TABLE(row,3));
end
fprintf('%s\n',repmat('-',1,59));
fprintf('%s%s       %3.0f\n\n\n','Total',repmat(' ',1,25),cols);