function [brainarea] = arearead_TT(session, cellname)
%

filename = ['TT_',session]; eval([filename])
monkey = parsession.monkey;
electrode = str2num(cellname(1));
locs = parsession.locations;
uprobe = parsession.uprobe;

if ~isempty(locs)
    if ~uprobe    
        ind = find(locs(:,1)==electrode);
        loc = locs(ind,:);
    else
        loc = locs;
    end
else
	loc=[];
end


if loc(4)==4
    brainarea = 'OFC'; 
elseif loc(4)==9
    brainarea = 'DLPFC';
elseif loc(4)==10
    brainarea = 'VLPFC';
else
    brainarea = 'other';
end



