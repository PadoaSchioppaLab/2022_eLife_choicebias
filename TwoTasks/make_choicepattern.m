function [pairlist, table01_all] = make_choicepattern(goodTrials)
% function [table01_all] = make_choicepattern(TrialRecord)
%
% goodTrials is a matrix with one row per trial and ngoods+3 columns, with
% [trialnumber, offers, chosen good, gotjuice] (only successful trials)
%
%

%
verbose = 0;
oneplot = 1;
disp_pooldir = 1;

%
% goodTrials = TrialRecord.goodTrials;
ngoods = size(goodTrials,2) - 3;

%
%compute choice patterns
offerlist = goodTrials(:,2:ngoods+1);
pairlist = [];
table01_all.bydir = {};
table01_all.pooldir = {};
for ig = 1:ngoods-1
	for jg = ig+1:ngoods
		ind = find(offerlist(:,ig) | offerlist(:,jg)); %FIXSB?
		if ~isempty(ind)
			pairlist = [pairlist; ig jg];

			%add forced choices
			ind = [ind; find((offerlist(:,ig) | offerlist(:,jg)) & sum(offerlist~=0,2)==1)];
			pairtrials = [offerlist(ind,[ig jg]), goodTrials(ind,end-1:end)];

			%
			%separate by left/right
			for isign = [-1 1]
				ind = find(sign(pairtrials(:,1))==isign | sign(pairtrials(:,2))==-isign);
				pairtrials_Aside = pairtrials(ind,:);
				%
				[offer, j, groups] = unique(pairtrials_Aside(:,1:2),'rows');
				choiz = pairtrials_Aside(:,3);
				%
				perc_B =[]; Ntrials = [];
				for i = 1:size(offer,1)
					ind = find(groups == i);
					perc_B(i,1) = length(find(choiz(ind) == jg))/length(ind);	%perc of ch b
					Ntrials(i,1) = length(ind);
				end
				table01 = [offer, perc_B, Ntrials];
				%
				%sort choices
				eps = 0.001;
				aux = abs(table01) + eps;
				[ratios, jnd] = sort(aux(:,2)./aux(:,1));
				%					[ratios, jnd] = sort(abs(table01(:,2))./abs(table01(:,1)));
				table01_all.bydir{size(pairlist,1),isign/2+1.5} = table01(jnd,:);
			end %for isign

			%
			%pooling directions
			[offer, j, groups] = unique(abs(pairtrials(:,1:2)),'rows');
			choiz = pairtrials(:,3);
			%
			perc_B =[]; Ntrials = [];
			for i = 1:size(offer,1)
				ind = find(groups == i);
				perc_B(i,1) = length(find(choiz(ind) == jg))/length(ind);	%perc of ch b
				Ntrials(i,1) = length(ind);
			end
			table01 = [offer, perc_B, Ntrials];
			%
			%sort choices
			eps = 0.001;
			aux = abs(table01) + eps;
			[ratios, jnd] = sort(aux(:,2)./aux(:,1));
			table01_all.pooldir{size(pairlist,1),1} = table01(jnd,:);

		end %if ~isempty
	end %for jg
end %for ig

%
if ~verbose
	return
end

npairs = size(table01_all.bydir,1);

%
%plot choice patterns
if (~oneplot)	figure; set(gcf,'position',[260 660 1000 290]), hold on, box on
else			figure; set(gcf,'position',[550 480 600 480]), hold on, box on
end

if ~oneplot
	for ipair = 1:npairs
		hs(ipair) = subplot(1,npairs,ipair); hold on, box on

		%make sure that same offer for both Asigns
		alloffers = [table01_all.bydir{ipair,1}(:,1:2); table01_all.bydir{ipair,2}(:,1:2)];
		alloffers = unique(abs(alloffers),'rows');
		[junk, jnd] = sort(alloffers(:,2)./alloffers(:,1));
		alloffers = alloffers(jnd,:);

		for isign = 1:2
			table01 = table01_all.bydir{ipair,isign};
			abs_off = abs(table01(:,1:2));
			[junk ind jnd] = intersect(abs_off, alloffers, 'rows');
			xx = jnd;
			yy = table01(ind,3);
			[junk ii] = sort(xx);
			xx = xx(ii);
			yy = yy(ii);
			hp1 = plot(xx, yy, 'o','markersize',8);
			hp2 = plot(xx, yy, '-');
			if		isign == 1
				set(hp1,'color',[0 .8 0],'markerfacecolor',[0 .8 0])
				set(hp2,'color',[0 .8 0])
			elseif	isign == 2
				set(hp1,'color',[1 0 0],	'markerfacecolor',[1 0 0])
				set(hp2,'color',[1 0 0])
			end
		end

		%cosmetics
		noffers = size(alloffers,1);
		axis([.5, noffers+.5, 0, 1]);
		for ioff = 1:noffers
			xticlab{ioff} = [num2str(alloffers(ioff,1)),':',num2str(alloffers(ioff,2))];
		end
		set(gca, 'xtick',[1:noffers], 'xticklabel',xticlab, 'fontsize',7)

	end

elseif oneplot
	if ~disp_pooldir
		%get alloffers
		alloffers = [];
		for ipair = 1:npairs
			for isign = 1:2
				alloffers = [alloffers; table01_all.bydir{ipair,isign}(:,1:2)];
			end
		end
		alloffers = unique(abs(alloffers),'rows');
		%
		%sort offers
		eps = 0.001;
		aux = alloffers + eps;
		[ratios, jnd] = sort(aux(:,2)./aux(:,1));
		alloffers = alloffers(jnd,:);

		%
		for ipair = 1:npairs
			for isign = 1:2
				table01 = table01_all.bydir{ipair,isign};
				abs_off = abs(table01(:,1:2));
				[junk ind jnd] = intersect(abs_off, alloffers, 'rows');
				xx = jnd;
				yy = table01(ind,3);
				[junk ii] = sort(xx);
				xx = xx(ii);
				yy = yy(ii);
				hp1 = plot(xx, yy, 'o','markersize',8);
				hp2 = plot(xx, yy, '-');
				if		isign == 1
					clr = [0 1 0]*ipair/npairs;
				elseif	isign == 2
					clr = [1 0 0]*ipair/npairs;
				end
				set(hp1, 'color',clr, 'markerfacecolor',clr)
				set(hp2, 'color',clr)
			end
		end

	elseif disp_pooldir
		%get alloffers
		alloffers = [];
		for ipair = 1:npairs
			alloffers = [alloffers; table01_all.pooldir{ipair,1}(:,1:2)];
		end
		alloffers = unique(abs(alloffers),'rows');
		%sort offers
		eps = 0.001;
		aux = alloffers + eps;
		[ratios, jnd] = sort(aux(:,2)./aux(:,1));
		alloffers = alloffers(jnd,:);

		%
		for ipair = 1:npairs
			table01 = table01_all.pooldir{ipair,1};
			abs_off = abs(table01(:,1:2));
			[junk ind jnd] = intersect(abs_off, alloffers, 'rows');
			xx = jnd;
			yy = table01(ind,3);
			[junk ii] = sort(xx);
			xx = xx(ii);
			yy = yy(ii);
			%
			%plot
			clr = [1 1 1]*(ipair-1)/(npairs);
			hp2 = plot(xx, yy, '-', 'linewidth',2, 'color',clr);
			if (npairs<5)
				plot(xx, yy, 'o', 'markersize',14-3*ipair, 'linewidth',2, 'color',clr, 'markerfacecolor',clr);
			else
				plot(xx, yy, 'o', 'markersize',16-2*ipair, 'linewidth',2, 'color',clr, 'markerfacecolor',clr);
			end
		end
	end

	%cosmetics
	noffers = size(alloffers,1);
	axis([.5, noffers+.5, 0, 1]);
	for ioff = 1:noffers
		xticlab{ioff} = [num2str(alloffers(ioff,1)),':',num2str(alloffers(ioff,2))];
	end
	set(gca, 'xtick',[1:noffers], 'xticklabel',xticlab, 'fontsize',9)
end %if oneplot
% toc