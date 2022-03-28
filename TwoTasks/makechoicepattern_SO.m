function [pairlist, table01_all] = makechoicepattern_SO(goodTrials)
% function [table01_all] = make_choicepattern(TrialRecord)
%
% goodTrials is a matrix with one row per trial and ngoods+3 columns, with
% [trialnumber, offers, chosen good, gotjuice] (only successful trials)
%
%

%

%
% goodTrials = TrialRecord.goodTrials;
ngoods = size(goodTrials,2) - 3;

%make goodTrials in the original format for uppper/lower hemifield plot
juiceBup = ones(size(goodTrials,1),1);
ind = goodTrials(:,ngoods+2)<=4;
juiceBup(ind) = -1;
goodTrials(:,end-1) = goodTrials(:,end); 
goodTrials(:,end) = juiceBup;

%
%compute choice patterns
offerlist = goodTrials(:,2:ngoods+1);
pairlist = [];
table01_all.byorder = {};
table01_all.poolord = {};
for ig = 1:ngoods-1
	for jg = ig+1:ngoods
		ind = find(offerlist(:,ig) & offerlist(:,jg));
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
				table01_all.byorder{size(pairlist,1),isign/2+1.5} = table01(jnd,:);
			end %for isign

			%
			%separate by upper/lower tgt hemifield
			for ihmf = [-1 1]
				ind = sign(pairtrials(:,end)) == -ihmf;
				pairtrials_Ahmf = pairtrials(ind,:);
				pairtrials_Ahmf(:,1:2) = abs(pairtrials_Ahmf(:,1:2));
				%
				[offer, j, groups] = unique(pairtrials_Ahmf(:,[1 2 4]),'rows');
				choiz = pairtrials_Ahmf(:,3);
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
				table01_all.bymov{size(pairlist,1),ihmf/2+1.5} = table01(jnd,:);
			end %for ihmf


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
			table01_all.poolord{size(pairlist,1),1} = table01(jnd,:);

		end %if ~isempty
	end %for jg
end %for ig
