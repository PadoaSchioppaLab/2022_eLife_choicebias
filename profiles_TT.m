function [profile] = profiles_TT(cellname)
%
% makes and saves the variable profile. profile is a structure with
% fileds:	.type (for which sides are collapsed)
%			.side (for which pairs and types are collapsed)
% the firing profile is computed for various allignments.
% this function calls profilemake_RM to compute any single
% profile.
%


% author: camillo, january 2004. revisions:
% may 2007: added profiles by offertype, trialtype
% november 2008: modified for RC
% december 2009: modified for DT
% december 2010: 'downgraded' to JC
% december 2016: modified for SO by SB
% october  2018: modified for TT by WS

% clear all
% tic
% cellname = 'V101203a11';
disp(['   ... computing profiles for cell ',cellname])

%keyboard

% mode = 'binned';
mode = 'smooth';
atleast_nntrials = 3;	%at least this many trials to include trialtype in analysis

%sides
orders = [-1,1]; % for SO
sides = [-1,+1]; % for JC
side_names = {'left','right'}; % for JC

session = cellname(1:8); readsession_TT %#ok<NASGU>
allignments = {'offeron', 'sacctgton', 'choicemade'};
% allignments = {'offeron', 'fixoff', 'outcome'};
% keyboard
allflags	= [flags.offeron, flags.sacctgton flags.choicemade ];
% allflags	= [flags.offeron, flags.fixoff flags.choicemade ];
othflags = [flags.fixon, flags.offeron, flags.offeroff, flags.offer2on, flags.offer2off,...
	flags.sacctgton, flags.choicemade , flags.outcome];

%load cell data, bhvParams
filename = [dirroot, cellname, '_data']; eval(['load ',filename])
filename = [dirroot, cellname, '_bhvParams']; eval(['load ',filename])
ngoods = size(sessionParams.goods,2);
if isequal(pairs_JC, pairs_SO)
    pairs = pairs_JC;
else
    disp('error: pairs are different in JC and SO');
end

npairs = size(pairs,1);

profile.specs.mode	= mode;
profile.specs.pairs	= pairs;
profile.specs.trial	= sessionParams.trial;

%in psyphydata convert possible_outcomes to flag_outcome
ind = ismember(psyphydata(:,3), flags.poss_outcomes);
psyphydata(ind,3) = flags.outcome;




for ipair = 1:npairs
    
    %
    % JC
	table01 = table01_all_JC.pooldir{ipair};
	hitdata = get_hitdata_TT(goodTrials_JC, pairlist_JC(ipair,:), ngoods, 'JC');

    %profile.*.all
	hits = hitdata(:,1);
	for ialign = 1:length(allignments)
		if		isequal(mode,'binned')
			prf = profilemake_binned(celldata, psyphydata, hits, flags, allflags(ialign), 10);
		elseif	isequal(mode,'smooth')
			prf = profilemake_smooth_TT(celldata, psyphydata, hits, flags, allflags(ialign),'JC');
		end
		eval(['profile.JC.',pairs{ipair},'.all.',allignments{ialign},' = prf;'])
		
		% flagstimes
		[flagtimes]=get_flagtimes(psyphydata,hits, allflags(ialign),othflags);
		eval(['profile.JC.',pairs{ipair},'.all.flagtimes_',allignments{ialign},' = flagtimes;'])
	end
	
	
	%profile.*.offertype
	for itype = 1 : size(table01,1)
		ind = hitdata(:,2)==table01(itype,1) & hitdata(:,3)==table01(itype,2);
		hits_type = hitdata(ind,1);
		%
		for ialign = 1:length(allignments)
			if		isequal(mode,'binned')
				prf = profilemake_binned(celldata, psyphydata, hits_type, flags, allflags(ialign), 20);
			elseif	isequal(mode,'smooth')
				%keyboard
				prf = profilemake_smooth_TT(celldata, psyphydata, hits_type, flags, allflags(ialign),'JC');
			end
			
			if (table01(itype,1)<10), auxstr1 = [strcat('0'),num2str(table01(itype,1))]; else auxstr1 = num2str(table01(itype,1)); end
			if (table01(itype,2)<10), auxstr2 = [strcat('0'),num2str(table01(itype,2))]; else auxstr2 = num2str(table01(itype,2)); end
			type_name = [pairs{ipair}(1), auxstr1, pairs{ipair}(2), auxstr2];
			
			eval(['profile.JC.',pairs{ipair},'.offertype.',type_name,'.',allignments{ialign},' = prf;'])
			
			% flagstimes
			[flagtimes]=get_flagtimes(psyphydata,hits_type, allflags(ialign),othflags);
			eval(['profile.JC.',pairs{ipair},'.offertype.',type_name,'.flagtimes_',allignments{ialign},' = flagtimes;'])			
		end
    end
    
    %profile.*.trialtype
    hitchoice = hitdata(:,4).*hitdata(:,5);
    for itype = 1 : size(table01,1)
        for ichoice = [1,-1]  % 1 for A in pair (AB)
            ind = hitdata(:,2)==table01(itype,1) & hitdata(:,3)==table01(itype,2) & hitchoice==ichoice;
            hits_type = hitdata(ind,1);
            if size(hits_type,1)>atleast_nntrials
                for ialign = 1:length(allignments)
                    if		isequal(mode,'binned')
                        prf = profilemake_binned(celldata, psyphydata, hits_type, flags, allflags(ialign), 10);
                    elseif	isequal(mode,'smooth')
                        prf = profilemake_smooth_TT(celldata, psyphydata, hits_type, flags, allflags(ialign),'JC');
                    end
                    if (table01(itype,1)<10), auxstr1 = [strcat('0'),num2str(table01(itype,1))]; else auxstr1 = num2str(table01(itype,1)); end
                    if (table01(itype,2)<10), auxstr2 = [strcat('0'),num2str(table01(itype,2))]; else auxstr2 = num2str(table01(itype,2)); end
                    if (ichoice==1), auxstr3 = pairs{ipair}(1); elseif (ichoice==-1), auxstr3 = pairs{ipair}(2); else dummy; end
                    type_name = [pairs{ipair}(1), auxstr1, pairs{ipair}(2), auxstr2, '_', auxstr3];
                    eval(['profile.JC.',pairs{ipair},'.trialtype.',type_name,'.',allignments{ialign},' = prf;'])       
                    % flagstimes
                    [flagtimes]=get_flagtimes(psyphydata,hits_type, allflags(ialign),othflags);
                    eval(['profile.JC.',pairs{ipair},'.trialtype.',type_name,'.flagtimes_',allignments{ialign},' = flagtimes;'])
                end
            end %if atleast_nntrials
        end %for ichoice
    end

    %profile.*.side
    crits = {'juice','choice'};
    critinds = [4,5];
    for icrit = 1:2
        crit = crits{icrit};
        critind = critinds(icrit);
        for iside = 1:2
            side = sides(iside);
            ind = hitdata(:,critind)==side;
            hitdata_side = hitdata(ind,:);
            hits_side = hitdata_side(:,1);
            %
            for ialign = 1:length(allignments)
                if		isequal(mode,'binned')
                    prf = profilemake_binned(celldat, psyphydata, hits_side, flags, allflags(ialign), 10);
                elseif	isequal(mode,'smooth')
                	prf = profilemake_smooth_TT(celldata, psyphydata, hits_side, flags, allflags(ialign),'JC');
                end
                eval(['profile.JC.',pairs{ipair},'.side.',crit,'.',side_names{iside},'.',allignments{ialign},' = prf;'])
                % flagstimes
                [flagtimes]=get_flagtimes(psyphydata,hits_type, allflags(ialign),othflags);
                eval(['profile.JC.',pairs{ipair},'.side.',crit,'.',side_names{iside},'.flagtimes_',allignments{ialign},' = flagtimes;'])
            end
        end	%for iside
    end %for icrit
       
    
    %
    % SO
	table01 = table01_all_SO.pooldir{ipair};
	hitdata = get_hitdata_TT(goodTrials_SO, pairlist_SO(ipair,:), ngoods, 'SO');
	
	%profile.*.all
	hits = hitdata(:,1);
	for ialign = 1:length(allignments)
		if		isequal(mode,'binned')
			prf = profilemake_binned(celldata, psyphydata, hits, flags, allflags(ialign), 10);
		elseif	isequal(mode,'smooth')
			prf = profilemake_smooth_TT(celldata, psyphydata, hits, flags, allflags(ialign),'SO');
		end
		eval(['profile.SO.',pairs{ipair},'.all.',allignments{ialign},' = prf;'])
		
		% flagstimes
		[flagtimes]=get_flagtimes(psyphydata,hits, allflags(ialign),othflags);
		eval(['profile.SO.',pairs{ipair},'.all.flagtimes_',allignments{ialign},' = flagtimes;'])
	end
	
	
	%profile.*.offertype
	for itype = 1 : size(table01,1)
		ind = hitdata(:,2)==table01(itype,1) & hitdata(:,3)==table01(itype,2);
		hits_type = hitdata(ind,1);
		%
		for ialign = 1:length(allignments)
			if		isequal(mode,'binned')
				prf = profilemake_binned(celldata, psyphydata, hits_type, flags, allflags(ialign), 20);
			elseif	isequal(mode,'smooth')
				%keyboard
				prf = profilemake_smooth_TT(celldata, psyphydata, hits_type, flags, allflags(ialign),'SO');
			end
			
			if (table01(itype,1)<10), auxstr1 = [strcat('0'),num2str(table01(itype,1))]; else auxstr1 = num2str(table01(itype,1)); end
			if (table01(itype,2)<10), auxstr2 = [strcat('0'),num2str(table01(itype,2))]; else auxstr2 = num2str(table01(itype,2)); end
			type_name = [pairs{ipair}(1), auxstr1, pairs{ipair}(2), auxstr2];
			
			eval(['profile.SO.',pairs{ipair},'.offertype.',type_name,'.',allignments{ialign},' = prf;'])
			
			% flagstimes
			[flagtimes]=get_flagtimes(psyphydata,hits_type, allflags(ialign),othflags);
			eval(['profile.SO.',pairs{ipair},'.offertype.',type_name,'.flagtimes_',allignments{ialign},' = flagtimes;'])
			
		end
	end
	
	
	for criter=1:6
			%	keyboard
		%FILTER FOR PROFILE
		%1: Choice by ChoosenID: 5: 1=Choosen A OR -1=Choosen B
		%2: Choice by Order: 5*6: 1=Choosen First OR -1=Choosen Second
		%3: Choice by Side: 5*4: 1=Choosen Left OR -1=Choosen Right
		%4: Offer by Order: 6: 1=Offer_AB OR -1=Offer_BA

		switch criter
			case 1
				hitchoice = hitdata(:,5);
				order_names = {'ChosenJuice.B','ChosenJuice.A'};
			case 2
				hitchoice = hitdata(:,5).*hitdata(:,6);
				order_names = {'ChosenOrder.Second','ChosenOrder.First'};
			case 3
				hitchoice = hitdata(:,5).*hitdata(:,4);
				order_names = {'ChosenSide.Right','ChosenSide.Left'};	
			case 4
				hitchoice = hitdata(:,6);
				order_names = {'OfferBA.all','OfferAB.all'};
			case 5
				for n=1:size(hitdata,1)
					if hitdata(n,6)==-1
					hitchoice(n,:)=hitdata(n,5);
					else
						hitchoice(n,:)=NaN;
					end	
				end
				order_names = {'OfferBA.ChB','OfferBA.ChA'};
			case 6
				for n=1:size(hitdata,1)
					if hitdata(n,6)==1
					hitchoice(n,:)=hitdata(n,5);
					else
						hitchoice(n,:)=NaN;
					end	
				end
				order_names = {'OfferAB.ChB','OfferAB.ChA'};	
		end
		
		
		
		% SAVE SPECS  for profile plot % SB
		profile.specs.allignments=allignments;
		eval(['profile.specs.criteria{', num2str(criter) ,'}','=order_names;']);
		
		%profile.*.ORDER and *.trialtype
		for iorder= 1:2
			mask=hitchoice==orders(iorder);
			hit_order = 	hitdata(mask,:);
			for ialign = 1:length(allignments)
				%  BY ORDER
				if		isequal(mode,'binned')
					prf = profilemake_binned(celldata, psyphydata, hit_order, flags, allflags(ialign), 10);
				elseif	isequal(mode,'smooth')
					prf = profilemake_smooth_TT(celldata, psyphydata, hit_order, flags, allflags(ialign),'SO');
				end
				eval(['profile.SO.',pairs{ipair},'.',order_names{iorder},'.',allignments{ialign},' = prf;'])
				
				% flagstimes
				[flagtimes]=get_flagtimes(psyphydata,hit_order, allflags(ialign),othflags);
				eval(['profile.SO.',pairs{ipair},'.',order_names{iorder},'.flagtimes_',allignments{ialign},' = flagtimes;'])
				
				% BY TRIAL TYPE
				for itype = 1 : size(table01,1)
					ind = hit_order(:,2)==table01(itype,1) & hit_order(:,3)==table01(itype,2);
					if sum(ind)>=atleast_nntrials
						hits_type = hit_order(ind,1);
						if		isequal(mode,'binned')
							prf = profilemake_binned(celldata, psyphydata, hits_type, flags, allflags(ialign), 10);
						elseif	isequal(mode,'smooth')
							prf = profilemake_smooth_TT(celldata, psyphydata, hits_type, flags, allflags(ialign),'SO');
						end
						
						if (table01(itype,1)<10), auxstr1 = [strcat('0'),num2str(table01(itype,1))]; else auxstr1 = num2str(table01(itype,1)); end
						if (table01(itype,2)<10), auxstr2 = [strcat('0'),num2str(table01(itype,2))]; else auxstr2 = num2str(table01(itype,2)); end
						type_name = [pairs{ipair}(1), auxstr1, pairs{ipair}(2), auxstr2];
						eval(['profile.SO.',pairs{ipair},'.',order_names{iorder},'.',type_name,'.',allignments{ialign},' = prf;'])

						% flagstimes
						[flagtimes]=get_flagtimes(psyphydata,hits_type, allflags(ialign),othflags);
						eval(['profile.SO.',pairs{ipair},'.',order_names{iorder},'.',type_name,'.flagtimes_',allignments{ialign},' = flagtimes;'])
					end
					% trial type in specs
					eval(['profile.specs.trialtypes.',type_name,' = [];'])						
                end				
			end	 %for ialign			
		end %for iorder		
	end %end criter
	
end %for ipair

%save profile
filename = [dirroot,cellname,'_profiles'];
eval(['save ',filename,' profile'])



%%%%%%%  functions  %%%%%%%

function [profile] = profilemake_binned(celldata, psyphydata, trials, flags, allflag, pace)
% computes the spikeing profile (in spikes per second) in a "binned" way.
% (pace is an input)
if		(allflag==flags.offeron),	llim = -750;  rlim =3500;	%offer on
elseif	(allflag==flags.choicemade),	llim =-500;  rlim = 2500;	%choicemade 
	%elseif	(allflag==flags.outcome),	llim =-1000;  rlim =2000;	%juice del
end
%
if isempty(trials)
	disp('no trials here')
	a = [llim:pace:rlim-pace]';
	profile = [a,zeros(size(a,1),1)];
	return
end
%
ind = ismember(celldata(:,2),trials);
sts = celldata(ind,:);						% sts is for spike time stamps
ind = ismember(psyphydata(:,2),trials);
pts = psyphydata(ind,:);					% pts is for psypchphysics time stamps
ntrials = length(trials);
allspikes = [];
%
for ii = 1:ntrials
	currtrial = trials(ii);
	try
		ssts = sts(sts(:,2)==currtrial, :);
		ppts = pts(pts(:,2)==currtrial, :);
		%reallign with allflag
		alligntime = ppts(ppts(:,3)==allflag,1);
		ssts(:,1) = ssts(:,1) - alligntime;
		
		ind2 = ssts(:,1)>llim & ssts(:,1)<rlim;
		allspikes = [allspikes; ssts(ind2,1)];
	catch
		disp(['some problem at trial number ',num2str(currtrial)])
	end
end
%
profile = [];
for ileft = llim : pace : rlim-pace
	iright = ileft + pace;
	nsp = length(find(allspikes>ileft & allspikes<iright));
	profile = [profile;(iright+ileft)/2,nsp];
end
profile(:,2) = profile(:,2)/ntrials/pace*1000;

