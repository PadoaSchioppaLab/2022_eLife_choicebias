function [psyphycell_newforTT] = sigmoidfit_TT_OrdChHyst(cellname,fitt,HystType,spaceType,varargin)
% this function does sigmoidal fitting with consideration of order and
% choice hysteresis
% this function is called by pop_behav_ana_summary_TT.m
%
% author: camillo, september 2004.
% revisions:
% june 2006: sigmoid is computed cell by cell + common value scale for JC3
% march 2007: added startpoint to fit sigmoid function
% november 2008: adapted for multiple choice patterns (RC). removed opportunity value
% october 2018: adapted for two tasks (JC and SO) -WS
% september 2019: add variable of hysteresis in the fitting for two tasks -WS

if isempty(varargin) verbose = 0;
else verbose = varargin{1};
end

if ~verbose
% 	clearvars
% 	tic
% 	cellname = 'G180912c11';
%     fitt = 'probit';
%     HystType = 'Only'; % 'Only': SO on SO or JC on JC; 'Both': SO and JC on either SO or JC
% 	verbose = 1;
end

disp(['   ... fitting sigmoid (by cell): ',cellname]);

fl=0; leg{1}='AB';leg{2}='BA'; leg{3}='JC'; %; leg{2}='AB'; leg{3}='BA';

fl=1; %fl+1;

space = spaceType;      % 'log' or 'linear'
dx = .005;			% for Riemann's integral
epsilon = 0.0001;	% for integration domain (how small the distr. before we ignore)

%load data
session = cellname(1:8); readsession_TT
filename = [dirroot, cellname, '_data']; eval(['load ',filename])
filename = [dirroot, cellname, '_bhvParams']; eval(['load ',filename])
npairs_JC = size(pairs_JC,1);
npairs_SO = size(pairs_SO,1);
if npairs_JC~=npairs_SO
    disp('error: pair numbers are different between JC and SO');
else
    npairs = npairs_JC;
end
if pairlist_JC~=pairlist_SO
    disp('error: pair lists are different between JC and SO');
else
    pairlist = pairlist_JC;
end

%probability ratios
prob_ratios = zeros(npairs,1);
for ipair = 1:npairs
	prob_ratios(ipair) = sessionParams.goods(pairlist(ipair,1)).probability / sessionParams.goods(pairlist(ipair,2)).probability;
end

%alloffers_num
%alloffers_num_JC
alloffers_num_JC = [];
for ipair = 1:npairs
	table01_JC = table01_all_JC.pooldir{ipair};
	alloffers_num_JC = [alloffers_num_JC; table01_JC(:,1:2)];
end
alloffers_num_JC = unique(abs(alloffers_num_JC),'rows');
eps = 0.001;
aux = alloffers_num_JC + eps;
%keyboard
[~, jnd] = sort(aux(:,2)./aux(:,1));
alloffers_num_JC = alloffers_num_JC(jnd,:);
%
%alloffers_num_SO
alloffers_num_SO = [];
for ipair = 1:npairs
	table01_SO = table01_all_SO.pooldir{ipair};
	alloffers_num_SO = [alloffers_num_SO; table01_SO(:,1:2)];
end
alloffers_num_SO = unique(abs(alloffers_num_SO),'rows');
eps = 0.001;
aux = alloffers_num_SO + eps;
%keyboard
[~, jnd] = sort(aux(:,2)./aux(:,1));
alloffers_num_SO = alloffers_num_SO(jnd,:);
%
if alloffers_num_JC~=alloffers_num_SO
    disp('error: all offers are different between JC and SO');
else
    alloffers_num = alloffers_num_JC;
end



for ipair = 1:npairs
	for dirswitch=1
		
        %
        % SO only (SO hysteresis on SO)
        %
        % column in goodTrials_SO: trial#, offerA, offerB, posA, chosen juice
        goodTrials = goodTrials_SO;
        goodTrials(:,4) = goodTrials_SO(:,5);
        goodTrials(:,5) = goodTrials_SO(:,4);
        nhits = size(goodTrials,1); %#ok<NODEF>
        aux = goodTrials(:,2:3);
        aux2 = nan(nhits,1);	 % chosen order for SO % 1 for second, -1 for first
        for i=1:nhits
            aux2(i)=sign(aux(i,goodTrials(i,4))); 
        end
        % trial#, offerA, offerB, chosen order(SO)
        hitinfo = [goodTrials(:,1:3), aux2]; 
        ordA = sign(sign(hitinfo(:,2) - hitinfo(:,3))); % 1 for BA and -1 for AB
        % trial#, offerA, offerB, order(SO), chosen order(SO), chosen juice
        hitdata_SO = [abs(hitinfo(:,1:3)), ordA, hitinfo(:,4)]; 
        hitdata_SO(:,6) = -ones(nhits,1); % 1 for B and -1 for A
        hitdata_SO(ordA.*hitinfo(:,4)<0,end) = 1;
        %
        % add choice hysteresis
        indA = hitdata_SO(:,6)==-1;
        trA  = hitdata_SO(indA,1);
        indB = hitdata_SO(:,6)==1;
        trB  = hitdata_SO(indB,1);
        hitdata_SO(:,7) = 0; % choice of previous trial: -1 for A, 1 for B, 0 for else
        hitdata_SO(ismember(hitdata_SO(:,1),trA+1),7) = -1; 
        hitdata_SO(ismember(hitdata_SO(:,1),trB+1),7) = 1; 
        % 
        %
        yy_SO  = (hitdata_SO(:,6)+1)/2; % choice: 1 for B and 0 for A;
        if isequal(space,'linear')
            xx1_SO = [hitdata_SO(:,3) , hitdata_SO(:,2)]; % [#B , #A]
        elseif	isequal(space,'log')
            xx1_SO = hitdata_SO(:,3)./hitdata_SO(:,2); % #B/#A
			xx1_SO = log(xx1_SO); % log(#B/#A)
		end
        xx2_SO = hitdata_SO(:,4); % order: 1 for BA and -1 for AB
        xx3_SO = hitdata_SO(:,7); % choice hysteresis: 1 for B choice, -1 for A choice and 0 for others
        ind_forcedA_SO = logical(hitdata_SO(:,2) & ~hitdata_SO(:,3));
        ind_forcedB_SO = logical(~hitdata_SO(:,2) & hitdata_SO(:,3));
        ii_SO = ~ind_forcedA_SO & ~ind_forcedB_SO; % index of non-forced choice trials
        
        %
        % JC only (JC hysteresis on JC)
        %
		% column in goodTrials_JC: trial#, offerA, offerB, chosen juice, gotjuice
        goodTrials = goodTrials_JC; 
		nhits = size(goodTrials,1); %#ok<NODEF>
        aux = goodTrials(:,2:3);
        aux2 = nan(nhits,1);	 % chosen side for JC 
        for i=1:nhits
            aux2(i)=sign(aux(i,goodTrials(i,4))); 
        end
        % trial#, offerA, offerB, chosen side(JC) chosen juice
        hitdata_JC = [abs(goodTrials(:,1:3)), aux2, goodTrials(:,4)]; 
		%
        % add choice hysteresis
        indA = hitdata_JC(:,5)==1;
        trA  = hitdata_JC(indA,1);
        indB = hitdata_JC(:,5)==2;
        trB  = hitdata_JC(indB,1);
        hitdata_JC(:,6) = 0; % choice of previous trial: -1 for A, 1 for B, 0 for else
        hitdata_JC(ismember(hitdata_JC(:,1),trA+1),6) = -1; 
        hitdata_JC(ismember(hitdata_JC(:,1),trB+1),6) = 1; 
        %
        yy_JC  = hitdata_JC(:,5)-1; % choice: 1 for B and 0 for A;
        if isequal(space,'linear')
            xx1_JC = [hitdata_JC(:,3), hitdata_JC(:,2)]; % [#B, #A]
        elseif	isequal(space,'log')
            xx1_JC = hitdata_JC(:,3)./hitdata_JC(:,2); % #B/#A
			xx1_JC = log(xx1_JC); % log(#B/#A)
        end
        xx2_JC = zeros(size(xx1_JC,1),1); % order: 0 for JC
        xx3_JC = hitdata_JC(:,6); % choice hysteresis: 1 for B choice, -1 for A choice and 0 for others
        ind_forcedA_JC = logical(hitdata_JC(:,2) & ~hitdata_JC(:,3));
        ind_forcedB_JC = logical(~hitdata_JC(:,2) & hitdata_JC(:,3));
        ii_JC = ~ind_forcedA_JC & ~ind_forcedB_JC; % index of non-forced choice trials
        
        
        %
        % JC and SO both (both JC and SO hysteresis on either SO or JC)
        %
		%
        preChoice_both = [[hitdata_SO(:,1);hitdata_JC(:,1)],[(hitdata_SO(:,6)+1)/2;hitdata_JC(:,5)-1]];  % choice: 1 for B and 0 for A;
        indA = preChoice_both(:,2)==0;
        trA  = preChoice_both(indA,1);
        indB = preChoice_both(:,2)==1;
        trB  = preChoice_both(indB,1);
        %
        HystBoth_SO = zeros(size(hitdata_SO,1),1); % choice of previous trial(both JC and SO): -1 for A, 1 for B, 0 for else
        HystBoth_SO(ismember(hitdata_SO(:,1),trA+1),1) = -1; 
        HystBoth_SO(ismember(hitdata_SO(:,1),trB+1),1) = 1; 
        %
        HystBoth_JC = zeros(size(hitdata_JC,1),1); % choice of previous trial(both JC and SO): -1 for A, 1 for B, 0 for else
        HystBoth_JC(ismember(hitdata_JC(:,1),trA+1),1) = -1; 
        HystBoth_JC(ismember(hitdata_JC(:,1),trB+1),1) = 1; 
        
        %
        yy_both = [yy_SO;yy_JC]; % choice: 1 for B and 0 for A;
        xx1_both = [xx1_SO;xx1_JC]; % log(#B/#A)
        xx2_both = [xx2_SO;xx2_JC]; % order: 0 for JC, 1 for BA and -1 for AB
        xx3_both = [HystBoth_SO;HystBoth_JC]; % one step choice hysteresis: 1 for B choice, -1 for A choice and 0 for others
        xx4_both = [ones(size(hitdata_SO,1),1);zeros(size(hitdata_JC,1),1)]; % task type: 1 for SO and 0 for JC
        ii_both  = [ii_SO;ii_JC];
        
        % Fitting the behavioral curve
		% CHOOSING THE FITTING FUNCTION
		if strcmp(fitt,'normcdf')
            % NOT USING ANYMORE!! -WS 09/2019
            disp('fitt is normcdf which is not used anymore!');
            keyboard
        else
            %
            % SO
            %
			%keyboard
			try				
				% full fitting: choice ~ offertype + order + choice hysteresis
                if isequal(HystType,'Both');xx3_SO = HystBoth_SO; end
                if isequal(space,'linear')
                    [~ , devv, stats] = glmfit([xx1_SO(ii_SO,:),xx2_SO(ii_SO),xx3_SO(ii_SO)], yy_SO(ii_SO), 'binomial', 'link',fitt, 'constant','off');
                    psyphycell_newforTT.SO.OrdChHyst.sigmoidfit.vars = {'#B','#A','order','chHyst'};
                elseif	isequal(space,'log')
                    [~ , devv, stats] = glmfit([xx1_SO(ii_SO),xx2_SO(ii_SO),xx3_SO(ii_SO)], yy_SO(ii_SO), 'binomial', 'link',fitt, 'constant','on');
                    psyphycell_newforTT.SO.OrdChHyst.sigmoidfit.vars = {'constant','log(#B/#A)','order','chHyst'};
                end
                psyphycell_newforTT.SO.OrdChHyst.sigmoidfit.beta = stats.beta;
                psyphycell_newforTT.SO.OrdChHyst.sigmoidfit.se = stats.se;
                psyphycell_newforTT.SO.OrdChHyst.sigmoidfit.pval = stats.p;
                psyphycell_newforTT.SO.OrdChHyst.sigmoidfit.devv = devv;
                % tables
                ind_ABpreA = hitdata_SO(:,4)==-1 & hitdata_SO(:,7)==-1;
                pairtrials_SO = [hitdata_SO(ind_ABpreA,2),hitdata_SO(ind_ABpreA,3), yy_SO(ind_ABpreA,:)]; table01 = get_table01(pairtrials_SO);
                psyphycell_newforTT.SO.OrdChHyst.table01.ABpreA = table01;
                %
                ind_ABpreB = hitdata_SO(:,4)==-1 & hitdata_SO(:,7)==1;
                pairtrials_SO = [hitdata_SO(ind_ABpreB,2),hitdata_SO(ind_ABpreB,3), yy_SO(ind_ABpreB,:)]; table01 = get_table01(pairtrials_SO);
                psyphycell_newforTT.SO.OrdChHyst.table01.ABpreB = table01;
                %
                ind_ABpreO = hitdata_SO(:,4)==-1 & hitdata_SO(:,7)==0;
                pairtrials_SO = [hitdata_SO(ind_ABpreO,2),hitdata_SO(ind_ABpreO,3), yy_SO(ind_ABpreO,:)]; table01 = get_table01(pairtrials_SO);
                psyphycell_newforTT.SO.OrdChHyst.table01.ABpreO = table01;
                %
                ind_BApreA = hitdata_SO(:,4)==-1 & hitdata_SO(:,7)==-1;
                pairtrials_SO = [hitdata_SO(ind_BApreA,2),hitdata_SO(ind_BApreA,3), yy_SO(ind_BApreA,:)]; table01 = get_table01(pairtrials_SO);
                psyphycell_newforTT.SO.OrdChHyst.table01.BApreA = table01;
                %
                ind_BApreB = hitdata_SO(:,4)==-1 & hitdata_SO(:,7)==1;
                pairtrials_SO = [hitdata_SO(ind_BApreB,2),hitdata_SO(ind_BApreB,3), yy_SO(ind_BApreB,:)]; table01 = get_table01(pairtrials_SO);
                psyphycell_newforTT.SO.OrdChHyst.table01.BApreB = table01;
                %
                ind_BApreO = hitdata_SO(:,4)==-1 & hitdata_SO(:,7)==0;
                pairtrials_SO = [hitdata_SO(ind_BApreO,2),hitdata_SO(ind_BApreO,3), yy_SO(ind_BApreO,:)]; table01 = get_table01(pairtrials_SO);
                psyphycell_newforTT.SO.OrdChHyst.table01.BApreO = table01;
                %
                %
                % no hyst: choice ~ offertype + order
                if isequal(space,'linear')                   
                    [~ , devv, stats] = glmfit([xx1_SO(ii_SO,:),xx2_SO(ii_SO)], yy_SO(ii_SO), 'binomial', 'link',fitt, 'constant','off');
                    psyphycell_newforTT.SO.NonChHyst.sigmoidfit.vars = {'#B','#A','order'};
                elseif	isequal(space,'log')
                    [~ , devv, stats] = glmfit([xx1_SO(ii_SO),xx2_SO(ii_SO)], yy_SO(ii_SO), 'binomial', 'link',fitt, 'constant','on');
                    psyphycell_newforTT.SO.NonChHyst.sigmoidfit.vars = {'constant','log(#B/#A)','order'};
                end               
                psyphycell_newforTT.SO.NonChHyst.sigmoidfit.beta = stats.beta;
                psyphycell_newforTT.SO.NonChHyst.sigmoidfit.se = stats.se;
                psyphycell_newforTT.SO.NonChHyst.sigmoidfit.pval = stats.p;
                psyphycell_newforTT.SO.NonChHyst.sigmoidfit.devv = devv;
                
                % no hyst: choice ~ offertype (seperate fitting for AB and BA)
                ii_AB = (xx2_SO == -1) & ii_SO;
                ii_BA = (xx2_SO ==  1) & ii_SO;
                if isequal(space,'linear')                   
                    [~ , devv, stats] = glmfit([xx1_SO(ii_AB,:)], yy_SO(ii_AB), 'binomial', 'link',fitt, 'constant','off');
                    psyphycell_newforTT.SOAB.NonChHyst.sigmoidfit.vars = {'#B','#A'};
                elseif	isequal(space,'log')
                    [~ , devv, stats] = glmfit([xx1_SO(ii_AB)], yy_SO(ii_AB), 'binomial', 'link',fitt, 'constant','on');
                    psyphycell_newforTT.SOAB.NonChHyst.sigmoidfit.vars = {'constant','log(#B/#A)'};
                end                 
                psyphycell_newforTT.SOAB.NonChHyst.sigmoidfit.beta = stats.beta;
                psyphycell_newforTT.SOAB.NonChHyst.sigmoidfit.se = stats.se;
                psyphycell_newforTT.SOAB.NonChHyst.sigmoidfit.pval = stats.p;
                psyphycell_newforTT.SOAB.NonChHyst.sigmoidfit.devv = devv;
                if isequal(space,'linear')                   
                    [~ , devv, stats] = glmfit([xx1_SO(ii_BA,:)], yy_SO(ii_BA), 'binomial', 'link',fitt, 'constant','off');
                    psyphycell_newforTT.SOAB.NonChHyst.sigmoidfit.vars = {'#B','#A'};
                elseif	isequal(space,'log')
                    [~ , devv, stats] = glmfit([xx1_SO(ii_BA)], yy_SO(ii_BA), 'binomial', 'link',fitt, 'constant','on');
                    psyphycell_newforTT.SOAB.NonChHyst.sigmoidfit.vars = {'constant','log(#B/#A)'};
                end  
                psyphycell_newforTT.SOBA.NonChHyst.sigmoidfit.beta = stats.beta;
                psyphycell_newforTT.SOBA.NonChHyst.sigmoidfit.se = stats.se;
                psyphycell_newforTT.SOBA.NonChHyst.sigmoidfit.pval = stats.p;
                psyphycell_newforTT.SOBA.NonChHyst.sigmoidfit.devv = devv;
               
            catch
				keyboard
            end
            
			
            %
            % JC
            %
			%keyboard
			try				
				% full fitting: choice ~ offertype + choice hysteresis
                if isequal(HystType,'Both');xx3_JC = HystBoth_JC;end
                if isequal(space,'linear')                   
                    [~ , devv, stats] = glmfit([xx1_JC(ii_JC,:),xx3_JC(ii_JC)], yy_JC(ii_JC), 'binomial', 'link',fitt, 'constant','off');
                    psyphycell_newforTT.JC.OrdChHyst.sigmoidfit.vars = {'#B','#A','chHyst'};
                elseif	isequal(space,'log')
                    [~ , devv, stats] = glmfit([xx1_JC(ii_JC),xx3_JC(ii_JC)], yy_JC(ii_JC), 'binomial', 'link',fitt, 'constant','on');
                    psyphycell_newforTT.JC.OrdChHyst.sigmoidfit.vars = {'constant','log(#B/#A)','chHyst'};
                end                  
                psyphycell_newforTT.JC.OrdChHyst.sigmoidfit.beta = stats.beta;
                psyphycell_newforTT.JC.OrdChHyst.sigmoidfit.se = stats.se;
                psyphycell_newforTT.JC.OrdChHyst.sigmoidfit.pval = stats.p;
                psyphycell_newforTT.JC.OrdChHyst.sigmoidfit.devv = devv;
                % tables
                %
                ind_preA = hitdata_JC(:,6)==-1;
                pairtrials_JC = [hitdata_JC(ind_preA,2),hitdata_JC(ind_preA,3), yy_JC(ind_preA,:)]; table01 = get_table01(pairtrials_JC);
                psyphycell_newforTT.JC.OrdChHyst.table01.preA = table01;
                %
                ind_preB = hitdata_JC(:,6)==1;
                pairtrials_JC = [hitdata_JC(ind_preB,2),hitdata_JC(ind_preB,3), yy_JC(ind_preB,:)]; table01 = get_table01(pairtrials_JC);
                psyphycell_newforTT.JC.OrdChHyst.table01.preB = table01;
                %
                ind_preO = hitdata_JC(:,6)==0;
                pairtrials_JC = [hitdata_JC(ind_preO,2),hitdata_JC(ind_preO,3), yy_JC(ind_preO,:)]; table01 = get_table01(pairtrials_JC);
                psyphycell_newforTT.JC.OrdChHyst.table01.preO = table01;
                %
                %
                % no hyst: choice ~ offertype
                if isequal(space,'linear')                   
                    [~ , devv, stats] = glmfit([xx1_JC(ii_JC,:)], yy_JC(ii_JC), 'binomial', 'link',fitt, 'constant','off');
                    psyphycell_newforTT.JC.NonChHyst.sigmoidfit.vars = {'#B','#A'};
                elseif	isequal(space,'log')
                    [~ , devv, stats] = glmfit([xx1_JC(ii_JC)], yy_JC(ii_JC), 'binomial', 'link',fitt, 'constant','on');
                    psyphycell_newforTT.JC.NonChHyst.sigmoidfit.vars = {'constant','log(#B/#A)'};
                end 
                psyphycell_newforTT.JC.NonChHyst.sigmoidfit.beta = stats.beta;
                psyphycell_newforTT.JC.NonChHyst.sigmoidfit.se = stats.se;
                psyphycell_newforTT.JC.NonChHyst.sigmoidfit.pval = stats.p;
                psyphycell_newforTT.JC.NonChHyst.sigmoidfit.devv = devv;
            catch
				keyboard
            end
			
			
			%
            % BOTH (assume the same width)
            %
            try
                % full fitting: choice ~ offertype + order + choice hysteresis + task type
                if isequal(space,'linear')                   
                    [~ , devv, stats] = glmfit([xx1_both(ii_both,:),xx2_both(ii_both),xx3_both(ii_both),xx4_both(ii_both)], yy_both(ii_both), 'binomial', 'link',fitt, 'constant','off');
                    psyphycell_newforTT.both.OrdChHyst.sigmoidfit.vars = {'#B','#A','order','chHyst','tasktype'};
                elseif	isequal(space,'log')
                    [~ , devv, stats] = glmfit([xx1_both(ii_both),xx2_both(ii_both),xx3_both(ii_both),xx4_both(ii_both)], yy_both(ii_both), 'binomial', 'link',fitt, 'constant','on');
                    psyphycell_newforTT.both.OrdChHyst.sigmoidfit.vars = {'constant','log(#B/#A)','order','chHyst','tasktype'};
                end               
                psyphycell_newforTT.both.OrdChHyst.sigmoidfit.beta = stats.beta;
                psyphycell_newforTT.both.OrdChHyst.sigmoidfit.se = stats.se;
                psyphycell_newforTT.both.OrdChHyst.sigmoidfit.pval = stats.p;
                psyphycell_newforTT.both.OrdChHyst.sigmoidfit.devv = devv;
                
                % no hyst: choice ~ offertype + order + task type
                if isequal(space,'linear')                   
                    [~ , devv, stats] = glmfit([xx1_both(ii_both,:),xx2_both(ii_both),xx4_both(ii_both)], yy_both(ii_both), 'binomial', 'link',fitt, 'constant','off');
                    psyphycell_newforTT.both.NonChHyst.sigmoidfit.vars = {'#B','#A','order','tasktype'};
                elseif	isequal(space,'log')
                    [~ , devv, stats] = glmfit([xx1_both(ii_both),xx2_both(ii_both),xx4_both(ii_both)], yy_both(ii_both), 'binomial', 'link',fitt, 'constant','on');
                    psyphycell_newforTT.both.NonChHyst.sigmoidfit.vars = {'constant','log(#B/#A)','order','tasktype'};
                end
                psyphycell_newforTT.both.NonChHyst.sigmoidfit.beta = stats.beta;
                psyphycell_newforTT.both.NonChHyst.sigmoidfit.se = stats.se;
                psyphycell_newforTT.both.NonChHyst.sigmoidfit.pval = stats.p;
                psyphycell_newforTT.both.NonChHyst.sigmoidfit.devv = devv;
            catch
				keyboard
            end
            
        end		
	end % for dirswitch
end % for ipairs

end


%%
function [table01] = get_table01(pairtrials, dd)
%
[offer, ~, groups] = unique(abs(pairtrials(:,1:2)),'rows');
noffs = size(offer,1);
choiz = pairtrials(:,3);
%
perc_B = nan(noffs,1);
Ntrials = nan(noffs,1);
for i = 1:noffs
	ind = find(groups== i);
	perc_B(i,1) = length(find(choiz(ind)== 1))/length(ind);	%perc of ch b
	Ntrials(i,1) = length(ind);
end
table01 = [offer, perc_B, Ntrials];
%
%sort choices
eps = 0.001;
aux = abs(table01) + eps;
[~, jnd] = sort(aux(:,2)./aux(:,1));
table01 = table01(jnd,:);

if nargin==1, return, end
%include forced choice errors (FCE) in table01. currently works only for the analysis w/ all trials
for itype = 1:size(table01,1)
	if prod(table01(itype,1:2),2), continue, end
	ii = ismember(dd(:,[4,5]), table01(itype,[1,2]), 'rows');
	table01(itype,3) = sum(ii & dd(:,7)==1) / sum(ii & dd(:,7)<=2);
	%reverse if forced choice A (bcs of convention in dd(:,7))
	if ~table01(itype,2), table01(itype,3) = 1 - table01(itype,3); end
end
end