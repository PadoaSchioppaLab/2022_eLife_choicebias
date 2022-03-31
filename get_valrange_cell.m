function [valrange, valminmax_clas] = get_valrange_cell(parsession, clas, tuning, unival, JCorSO) %#ok<INUSL>
%
% computes the value range for offval cells and chval cells called by
%
% for taste cells, valminmax_clas and valrange are computed only in JC3
% (because in JC1 we dont know which juice is encoded). valminmax_class are
% the min and max offer values for the juice encoded by the cell, expressed
% in intrinsic units (uE)
%


%#ok<*NODEF>

[timewindows, ~] = get_timewindows_TT(JCorSO);

valminmax_clas = nan(1,2);	%we only compute for the encoded variable
valrange = nan;

if isequal(JCorSO,'JC')
    pairname = 'AB';
    tuning = tuning.JC;
elseif isequal(JCorSO,'SO')
    pairname = 'ABA';
    tuning = tuning.SO;
end
relvalue = unival;
%
%load typedata and remove trial types with <2% trials
eval(['typedata = tuning.',pairname,'.neuract.bytrialtype.',timewindows{1}{1},';'])
minperc = 2/100;
nt = typedata(:,end);
pt = nt/sum(nt);
ind = pt>minperc;
typedata = typedata(ind,:);

%offval cells
if ismember(clas,[1 2])
	valminmax = [min(typedata(:,1:2)); max(typedata(:,1:2))];
	valminmax(:,1) = valminmax(:,1)*relvalue;	%express all values in uB
	valminmax_clas = valminmax(:,clas);			%take only the clas
	valrange = diff(valminmax_clas);			%range as difference
	
	%chval cells
elseif clas==3
	offvals = typedata(:,1:2);
	chgoods = typedata(:,3);
	offvals(:,1) = offvals(:,1)*relvalue;		%express offvalA in uB
	%compute chvals
	ntypes = size(typedata,1);
	chvals = nan(ntypes,1);
	ind = find(chgoods== 1); chvals(ind,:) = offvals(ind,1);
	jnd = find(chgoods==-1); chvals(jnd,:) = offvals(jnd,2);
	valminmax_clas = [min(chvals(:,1)); max(chvals(:,1))];
	valrange = diff(valminmax_clas);			%range as difference
	
	%taste cells
elseif clas==4
	valminmax_clas = nan(1,2);
	valrange = nan;
end

