function BS = levlist2levlengths(GP,BS)
% Compute the levlengths info for a BasisSpec object
%
% Input
%   GP          a GraphPart object
%   BS          a BasisSpec object, without levlengths
%
% Output
%   BS          a BasisSpec object, with levlengths
%
%
%
% Copyright 2015 The Regents of the University of California
%
% Implemented by Jeff Irion (Adviser: Dr. Naoki Saito)



%% 0. Preliminaries

[levlist,~,c2f,description] = ExtractData(BS);

% allocate space
levlengths = zeros(length(levlist),1,class(GP.ind));


%% 1a. (COARSE-TO-FINE) Determine the length of each basis block specified by levlist
if c2f
    n = 0;
    for row = 1:length(levlist)
        % find the regionstart(s) of the next region
        IX = find( GP.rs(:,levlist(row))==n+1, 1, 'last' );

        levlengths(row) = GP.rs(IX+1,levlist(row))-GP.rs(IX,levlist(row));
        n = n+levlengths(row);
    end
    
    % get rid of blocks with length 0
    levlist(levlengths == 0) = [];
    levlengths(levlengths == 0) = [];


%% 1b. (FINE-TO-COARSE) Determine the length of each basis block specified by levlist
else
    n = 0;
    if isempty(GP.rsf2c)
        GP = FineToCoarse(GP);
    end
    
    for row = 1:length(levlist)
        IX = find( GP.rsf2c(:,levlist(row))==n+1, 1, 'last' );
        levlengths(row) = GP.rsf2c(IX+1,levlist(row))-GP.rsf2c(IX,levlist(row));
        n = n+levlengths(row);
    end
end


BS = BasisSpec(levlist,levlengths,c2f,description);


end
