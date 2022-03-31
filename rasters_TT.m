function [nplots] = rasters_TT(cellname)
%
% plots the rasters
%


disp(['   ... plotting rasters	of cell ',cellname])

session = cellname(1:8); readsession_TT %#ok<NASGU>

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
trialT = sessionParams.trial;
nplots=3;

allflag_SO = [flags.offeron, flags.choicemade];	sortflag = [-1 nan];	%-1 to sort trials by trialnumber
allflag_JC = [flags.offeron, flags.choicemade];	sortflag = [-1 nan];

%make raster plots, pair by pair
for ipair = 1:npairs
    %
    %JC
	table01 = table01_all_JC.pooldir{ipair};
	hitdata_JC = get_hitdata_TT(goodTrials_JC, pairlist_JC(ipair,:), ngoods,'JC');
	
	%type_names
	type_names = cell(size(table01,1),1);
	for itype = 1:size(table01,1)
		type_names{itype} = [num2str(table01(itype,1)), pairs{ipair}(1), ':', num2str(table01(itype,2)), pairs{ipair}(2)];
    end
	
    [~, hs] = rasterplot(celldata, psyphydata, hitdata_JC, table01, flags, allflag_JC, sortflag, type_names,'JC');
    %keyboard
    rasterplot_cosmetics(cellname, trialT, flags, allflag_JC, hs, 0); % order = 0; for JC
        
    %
    % SO
	table01 = table01_all_SO.pooldir{ipair};
	hitdata_SO = get_hitdata_TT(goodTrials_SO, pairlist_SO(ipair,:), ngoods,'SO');
	
	%type_names
	type_names = cell(size(table01,1),1);
	for itype = 1:size(table01,1)
		type_names{itype} = [num2str(table01(itype,1)), pairs{ipair}(1), ':', num2str(table01(itype,2)), pairs{ipair}(2)];
    end
    
	% RASTER PER OFFER ORDER
	for order=1:2
        ord=[-1 1];
        mask=(hitdata_SO(:,6)==ord(order));
        hitdataord=hitdata_SO(mask,:);
	
        [~, hs] = rasterplot(celldata, psyphydata, hitdataord, table01, flags, allflag_SO, sortflag, type_names,'SO');
        %keyboard
        
        rasterplot_cosmetics(cellname, trialT, flags, allflag_SO, hs,order);
    end
end
%toc



%%%%%%%%%%%%%%%%%%%%  FUNCTIONS  %%%%%%%%%%%%%%%%%%
function [hf, hs] = rasterplot(celldata, psyphydata, hitdata, table01, flags, allflag, sortflag, type_names, JCorSO)
psyphyonly = 0;
invertallign = 0;			%change for reverse the order of hits

nflags = length(allflag);
spikedata_align = cell(nflags,1);
psyphydata_align = cell(nflags,1);
hs = cell(nflags,1);

%
hf = figure; set(hf, 'Visible', 'on'); hold on;
for ialign = 1:nflags
	hs{ialign} = subplot(1,nflags,ialign); hold on;
end

sortref = ~isnan(sortflag);
sortflag = sortflag(~isnan(sortflag));

%in psyphydata convert possible_outcomes to flag_outcome
ind = ismember(psyphydata(:,3), flags.poss_outcomes);
psyphydata(ind,3) = flags.outcome;

	
%hits and spikedata
hits = hitdata(:,1);
nhits = size(hitdata,1);
ind = ismember(celldata(:,2),hits);
for ialign = 1:nflags
	spikedata_align{ialign} = celldata(ind,:);
	psyphydata_align{ialign} = psyphydata;
end


%reallign all psyphydata and spikedata
%try
for ialign = 1:nflags
	for ihit = 1:nhits
		trialnum = hits(ihit);
		allind = psyphydata(:,2)==trialnum & psyphydata(:,3)==allflag(ialign);
		alltime = psyphydata(allind,1);
		%alltime=min(alltime); %FIXSB
		hitind = find(psyphydata(:,2)==trialnum);
		if length(alltime)>1, [trialnum]; [alltime]; end
		psyphydata_align{ialign}(hitind,1) = psyphydata(hitind,1) - alltime;
		hitind = find(spikedata_align{ialign}(:,2)==trialnum);
		spikedata_align{ialign}(hitind,1) = spikedata_align{ialign}(hitind,1) - alltime;
	end
end
%catch
%	keyboard;
%end


%
y = -1.7;
yM = [];
yskipline = 4;
for itype = size(table01,1):-1:1
	y = y + yskipline;
	y0 = y;
	
	offtype = table01(itype,1:2);
	ind = hitdata(:,2)==offtype(1) & hitdata(:,3)==offtype(2);
	hitdata_type = hitdata(ind,1);
	hits_type = hitdata_type(:,1);
	
	%sort hits by sortflag
	hittable = [];
	psyphydata_ref = psyphydata_align{sortref};
	for ihit = 1:length(hits_type)
		trialnum = hits_type(ihit);
		if sortflag==-1		%trials sorted by trialnumber
			hittable(ihit,:) = [trialnum, trialnum];
		else
			flagind = psyphydata_ref(:,2)==trialnum & psyphydata_ref(:,3)==sortflag;
			flagtime = psyphydata_ref(flagind,1);
			if invertallign, flagtime = -flagtime; end
			hittable(ihit,:) = [trialnum, flagtime];
		end
	end
	%[hittable]
	[~,ind] = sort(hittable(:,2));
	sorthittable = hittable(ind,:);
	sorthittable = flipud(sorthittable);
	
	for ihit = 1:length(hits_type)
		%current trial, psyphy, spikes
		currtrial = sorthittable(ihit,1);
		%plot hits
		for ialign = 1:nflags
			subplot(hs{ialign})
			ind = psyphydata_align{ialign}(:,2)==currtrial;
			currpsyphy = psyphydata_align{ialign}(ind,:);
			ind = spikedata_align{ialign}(:,2)==currtrial;
			currspikes = spikedata_align{ialign}(ind,:);
			
			%psyphy events
			startTime		= currpsyphy(currpsyphy(:,3)==flags.trialstart, 1);
			trialendTime	= currpsyphy(currpsyphy(:,3)==flags.trialend, 1);
			
			%plot spikes
			if ialign==1
				y = y + .7;
			end
			if ~psyphyonly
 				ind2 = find(currspikes(:,1)>startTime & currspikes(:,1)<trialendTime+500);
				if ~isempty(ind2)
					lgtrial=trialendTime+500-startTime;
					toplot=zeros(1,lgtrial)*NaN;
					toplot(currspikes(ind2,1)-startTime)=currspikes(ind2,1);
					for n=1:4
					toplot(currspikes(ind2,1)-startTime+n)=currspikes(ind2,1)+n;
					end
					Y=ones(1,numel(toplot))*y;
					hp = plot(toplot,Y,'k','LineWidth',3);
% 					hp = plot(currspikes(ind2,1),y,'k.','markersize',3);	
				end
			end
			
			%plot psyphy events
			if ialign==1
				y = y + .3;
			end
			%
			ms = 2.5;
% 			ind3 = currpsyphy(:,3)==flags.trialstart;	%trial start
% 			plot(currpsyphy(ind3,1),y,'v','color',[.25 .25 .25],'markerfacecolor',[.25 .25 .25],'markersize',ms)
			
			ind3 = currpsyphy(:,3)==flags.fixon;		%fixon
			plot(currpsyphy(ind3,1),y,'>','color',[0.7 0.7 .7],'markerfacecolor',[0.7 0.7 .7],'markersize',ms)
			
% 			ind3 = currpsyphy(:,3)==flags.fixoff;		%fixof
% 			plot(currpsyphy(ind3,1),y,'<','color',[0.75 0.75 .75],'markerfacecolor',[0.75 0.75 .75],'markersize',ms)
% 			
			ind3 = currpsyphy(:,3)==flags.offeron;		%offer1 on
			plot(currpsyphy(ind3,1),y,'>','color',[0 0 .7],'markerfacecolor',[0 0 .7],'markersize',ms)
			
            if isequal(JCorSO,'SO')
                ind3 = currpsyphy(:,3)==flags.offeroff;		%offer off
                plot(currpsyphy(ind3,1),y,'<','color',[0 0 0.7],'markerfacecolor',[0 0 0.7],'markersize',ms)
            
                ind3 = currpsyphy(:,3)==flags.offer2on;		%offer2 on
                plot(currpsyphy(ind3,1),y,'>','color',[0 .7 0],'markerfacecolor',[0 .7 0],'markersize',ms)
			
                ind3 = currpsyphy(:,3)==flags.offer2off	;	%offer2 off
                plot(currpsyphy(ind3,1),y,'<','color',[0 .7 0],'markerfacecolor',[0 .7 0],'markersize',ms)
            end
            
			ind3 = currpsyphy(:,3)==flags.sacctgton;	%sacc tgts on
			plot(currpsyphy(ind3,1),y,'v','color',[.25 .25 .75],'markerfacecolor',[.25 .5 .75],'markersize',ms)
	        
            if isequal(JCorSO,'SO')
                ind3 = currpsyphy(:,3)==flags.fixoff;	%GO Signal
                plot(currpsyphy(ind3,1),y,'^','color',[.75 .5 .5],'markerfacecolor',[1 0 0],'markersize',ms)
            end
             
			ind3 = currpsyphy(:,3)==flags.choicemade;	%choice made
			plot(currpsyphy(ind3,1),y,'v','color',[1 0 1],'markerfacecolor',[1 0 1],'markersize',ms)
			
			ind3 = currpsyphy(:,3)==flags.outcome;	%trial outcome
			plot(currpsyphy(ind3,1),y,'v','color',[1 0 0],'markerfacecolor',[1 0 0],'markersize',ms)
			
			ind3 = currpsyphy(:,3)==flags.trialend;	%trial end
			plot(currpsyphy(ind3,1),y,'<','color',[.25 .25 .25],'markerfacecolor',[.25 .25 .25],'markersize',ms)
		
	
		
		end %for ialign
	end %for ihit
	
	yM(itype) = (y+y0)/2;
	
end %itype

for ialign = 1:nflags
	subplot(hs{ialign})
	set(gca,'ylim',[0 nhits + yskipline*size(table01,1)-1])
end

%add type_names
subplot(hs{1})
for itype = 1:size(table01,1)
	str = type_names{itype};
	ht = text(-900, yM(itype), str);
	% 	ht = text(-1900, yM(itype), str);
	set(ht,'rotation',90, 'horizontalalignment','center', 'fontname','arial', 'fontsize',8)
end



function [] = rasterplot_cosmetics(cellname, trialT, flags, allflag, hs, order)
do_highlights = 0;
% offeronT = trialT.offeron_time;
xrange = [];
if length(allflag)==1
	set(gca,'position',[0.15 0.06 0.80 0.865],'fontname','arial','fontsize', 8)
	if 		(allflag==flags.offeron)	set(gca,'xlim',[-1750 4500], 'xtick',[-1000:1000:4000])	%offeron
	elseif	(allflag==flags.outcome)	set(gca,'xlim',[-6000 1000], 'xtick',[-6000:1000:1000])	%juicedel
	end
else %if isequal(allflag,[flags.offeron flags.outcome]);
	subplot(hs{1}), set(gca,'xlim',[-1000 3500], 'xtick',[-1000:500:3500])			%offer1on
	xrange(1) = range(get(gca,'xlim'));
	subplot(hs{2}), set(gca,'xlim',[-1250 2750], 'xtick',[-1000:500:2500])			%choicemade
	xrange(2) = range(get(gca,'xlim'));
	
	% 	subplot(hs{3}), set(gca,'xlim',[-500 1000], 'xtick',[-500:500:1000])			%offer2on
	% 	xrange(3) = range(get(gca,'xlim'));
	% 	subplot(hs{3}), set(gca,'xlim',[-500 2000], 'xtick',[-500:500:2000])			%
	% 	xrange(3) = range(get(gca,'xlim'));
	%
	%keyboard
	%	pos1 = [0.08, .06, .8*xrange(1)/sum(xrange), .865];
	%	pos2 = [0.15+.8*xrange(1)/sum(xrange), .06, .8*xrange(2)/sum(xrange), .865];
	%	pos3 = [0.37+.8*xrange(1)/sum(xrange), .06, .8*xrange(2)/sum(xrange), .865];
	%	pos4 = [0.60+.8*xrange(1)/sum(xrange), .06, .8*xrange(2)/sum(xrange), .865];
	% 	subplot(hs{1}), set(gca, 'position',pos1, 'fontname','arial', 'fontsize', 8)
	% 	subplot(hs{2}), set(gca, 'position',pos2, 'fontname','arial', 'fontsize', 8)
	% 	subplot(hs{3}), set(gca, 'position',pos3, 'fontname','arial', 'fontsize', 8)
	% 	subplot(hs{4}), set(gca, 'position',pos4, 'fontname','arial', 'fontsize', 8)
	
end

for ialign = 1:length(allflag)
	subplot(hs{ialign})
	set(gca,'color',[0.95 0.95 0.95])
	grid on
	% 	set(gca,'ycolor',[.8 .8 .8])
	set(gca,'ycolor',ones(1,3))
end

%cellname
axes('position',[.7 .91 .5 .015])
text(0,3,['cell: ',cellname], 'fontname','arial', 'fontsize',9)
axis off

%position
%  set(gcf,'position',[25 10 520 960], 'PaperPositionMode','auto')
if order == 1
	title('AB')
set(gcf,'Units','normalized', 'position',[0 0 0.5 1], 'PaperPositionMode','auto')
elseif order == 2
	title('BA')
	set(gcf,'Units','normalized', 'position',[0.5 0 0.5 1], 'PaperPositionMode','auto')
elseif order == 0
	title('JC')
	set(gcf,'Units','normalized', 'position',[0.5 0 0.5 1], 'PaperPositionMode','auto')
end
%autoArrangeFigures();

%highlights
if do_highlights
	subplot(hs{1}),
	yyy=get(gca,'ylim');
	
	hp = fill([0 -500 -500 0],[0 0 yyy(2) yyy(2)],'w');
	set(hp, 'facealpha',.5, 'facecolor',[.7 .7 .7] ,'edgealpha',0)		%preoffer
	
	hp = fill([0 500 500 0],[0 0 yyy(2) yyy(2)],'w');
	set(hp, 'facealpha',.5, 'facecolor',[.3 .3 .7] ,'edgealpha',0)		%postoffer
	
	hp = fill([500 1000 1000 500],[0 0 yyy(2) yyy(2)],'w');
	set(hp, 'facealpha',.5, 'facecolor',[.7 .7 1] ,'edgealpha',0)		%latedelay
	
	subplot(hs{2})
	hp = fill([-1550 -1050 -1050 -1550],[0 0 yyy(2) yyy(2)],'w');
	set(hp, 'facealpha',.5, 'facecolor',[1 .3 .3] ,'edgealpha',0)		%prego
	
	hp = fill([-1050 -750 -750 -1050],[0 0 yyy(2) yyy(2)],'w');
	set(hp, 'facealpha',.5, 'facecolor',[.3 1 .3] ,'edgealpha',0)		%reactime
	
	hp = fill([0 -500 -500 0],[0 0 yyy(2) yyy(2)],'w');
	set(hp, 'facealpha',.5, 'facecolor',[.2 .75 .75] ,'edgealpha',0)	%prejuice
	
	hp = fill([0 500 500 0],[0 0 yyy(2) yyy(2)],'w');
	set(hp, 'facealpha',.5, 'facecolor',[.7 1 1] ,'edgealpha',0)		%postjuice
end

