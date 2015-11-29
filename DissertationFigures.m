function DissertationFigures
% Generate figures used in my dissertation
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



choice = menu('List of Reproducible Figures', 'Figure 2.8', 'Figure 2.9', 'Figure 3.2', ...
    'Figure 4.1', 'Figure 4.2', 'Figure 4.3', 'Figure 5.1', 'Figure 5.2', 'Figure 5.3', ...
    'Figure 5.4', 'Figure 5.5', 'Table 6.1', 'Figure 7.1', 'Figures 7.2, 7.3, & 7.4', ...
    'Figure 7.5', 'Figures 7.6, 7.7, & 7.8', 'Figure 8.1', 'Figures 8.2, 8.6a, and 8.7', ...
    'Figures 8.3, 8.6b, and 8.8', 'Figures 8.4, 8.6c, and 8.9', 'Figure 8.5', ...
    'Table 8.1 (takes a while)', 'Figure 9.2', 'Figure 9.3', 'Figure 9.4', 'Figure 9.5', ...
    'Figure 9.6', 'Figure 9.7', 'Figure 9.8', 'Figure 9.9', 'Table 10.1', 'Figure 10.1', ...
    'Figure 10.2', 'Figure 10.3', 'Figure 10.4', 'Figure 10.5', 'Figure 10.6', ...
    'Figure 10.7', 'Exit');

switch choice
    case 1
        JI_Figure2_8
        DissertationFigures
    case 2
        JI_Figure2_9
        DissertationFigures
    case 3
        JI_Figure3_2
        DissertationFigures
    case 4
        JI_Figure4_1
        DissertationFigures
    case 5
        JI_Figure4_2
        DissertationFigures
    case 6
        JI_Figure4_3
        DissertationFigures
    case 7
        JI_Figure5_1
        DissertationFigures
    case 8
        JI_Figure5_2
        DissertationFigures
    case 9
        JI_Figure5_3
        DissertationFigures
    case 10
        JI_Figure5_4
        DissertationFigures
    case 11
        JI_Figure5_5
        DissertationFigures
    case 12
        JI_Table6_1
        DissertationFigures
    case 13
        JI_Figure7_1
        DissertationFigures
    case 14
        JI_Figure7_2_3_4
        DissertationFigures
    case 15
        JI_Figure7_5
        DissertationFigures
    case 16
        JI_Figure7_6_7_8
        DissertationFigures
    case 17
        JI_Figure8_1
        DissertationFigures
    case 18
        JI_Figure8_2_6a_7
        DissertationFigures
    case 19
        JI_Figure8_3_6b_8
        DissertationFigures
    case 20
        JI_Figure8_4_6c_9
        DissertationFigures
    case 21
        JI_Figure8_5
        DissertationFigures
    case 22
        JI_Table8_1
        DissertationFigures
    case 23
        JI_Figure9_2
        DissertationFigures
    case 24
        JI_Figure9_3
        DissertationFigures
    case 25
        JI_Figure9_4
        DissertationFigures
    case 26
        JI_Figure9_5
        DissertationFigures
    case 27
        JI_Figure9_6
        DissertationFigures
    case 28
        JI_Figure9_7
        DissertationFigures
    case 29
        JI_Figure9_8
        DissertationFigures
    case 30
        JI_Figure9_9
        DissertationFigures
    case 31
        JI_Table10_1
        DissertationFigures
    case 32
        JI_Figure10_1
        DissertationFigures
    case 33
        JI_Figure10_2
        DissertationFigures
    case 34
        JI_Figure10_3
        DissertationFigures
    case 35
        JI_Figure10_4
        DissertationFigures
    case 36
        JI_Figure10_5
        DissertationFigures
    case 37
        JI_Figure10_6
        DissertationFigures
    case 38
        JI_Figure10_7
        DissertationFigures
    case 39
        close all
end


end