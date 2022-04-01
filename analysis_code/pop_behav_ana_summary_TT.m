% pop_behav_ana_summary_TT.m
%  
% this script summarizes the behavioral analysis of TT tasks of all sessions:
% 
% Relative value in JC and SO
% Steepness in JC and SO
% order bias in SO
% Choice hysteresis in JC and SO
% Reaction time in JC and SO
% Error rate in JC and SO
% ...

% Author:   WS - Nov 2019
% Revised:  WS - Feb 2020fit
% Revised:  WS - Jul 2020
% Revised:  WS - Oct 2020: find new examples
% Revised:  WS - Dec 2020: add analysis: "absolute" value; value range
% Revised:  WS - Apr 2021: add choice hysteresis in the fitting
% Revised:  WS - Feb 2022: add more option for bhv fitting: probit/logit; log ratio/linear diff

clearvars
close all

%parameters
brainarea = 'OFC';
% monkeys_ana = {'Juan'};  % 'Gervinho', 'Juan', 'both'
monkeys_ana = {'Juan', 'Gervinho'};  % 'Gervinho', 'Juan', 'both'
% monkeys_ana = {'both'};  % 'Gervinho', 'Juan', 'both'
nmonkeys = length(monkeys_ana);

doreacttime = 0;
doerrorrate = 0;
dovaluerange = 0;
dofractofA = 0;
domeanvalues = 0;

fittingwithCH  = 0;  % whether to do the probit regression with the term of choice hysteresis (both task types)

dologitorprobit = ''; % '_logit' or '_probit'

onlyplotpapers = 1;

plotjuicepairs = 0;

removesessions = 0; % remove based on dynamic range and saturation
remove_steepnessout = 1; % remove based on steepness outlier (1.5 IQR)

mintrialnum = 100; % 200, 160, 125, 100
atleast_nntrials = 2;

plotsessionexample = 0;
% plot session examples for each monkeys
examplecells = {'J190801b43','G181202a32'};
examplesessions = {'J190801b', 'G181202a'};
% examplecells = {'J190722b42', 'G181226a47'};
% examplesessions ={'J190722b', 'G181226a'};


% examplecells = {'G181202a32'};
% examplesessions = {'G181202a'};
if plotsessionexample  
    sigmoidfit_TT('G181226a47','probit',0,1);
    sigmoidfit_TT('J190722b42','probit',0,1);   
%     sigmoidfit_TT('G181202a32','probit',0,1);
%     sigmoidfit_TT('J190801b43','probit',0,1);  
%     % potential good ones
%     sigmoidfit_TT('G181016b11','probit',0,1);   
%     sigmoidfit_TT('G180911a42','probit',0,1);   
%     sigmoidfit_TT('G181023b41','probit',0,1);   
%     sigmoidfit_TT('G181127a42','probit',0,1);       
%     sigmoidfit_TT('G181227b43','probit',0,1);   
%     sigmoidfit_TT('G190113c33','probit',0,1); 
%     sigmoidfit_TT('G191225a31','probit',0,1); 

% good example to show preference bias and order bias
%     sigmoidfit_TT('G190103c41','probit',0,1); 
%     sigmoidfit_TT('G181226a47','probit',0,1); 
%     sigmoidfit_TT('G191212a12','probit',0,1); 
% 
%     sigmoidfit_TT('J190720a21','probit',0,1);
%     sigmoidfit_TT('J190730b42','probit',0,1);
%     sigmoidfit_TT('J190912a35','probit',0,1);
%     sigmoidfit_TT('J191108a43','probit',0,1); 
%     sigmoidfit_TT('J190719c42','probit',0,1);
%     sigmoidfit_TT('J190722a42','probit',0,1);
%     sigmoidfit_TT('J190730a22','probit',0,1);
%     sigmoidfit_TT('J190801b43','probit',0,1);
%     sigmoidfit_TT('J190723c42','probit',0,1);

% good example to show preference bias and order bias
%     sigmoidfit_TT('J190913a43','probit',0,1);
%     sigmoidfit_TT('J190805b33','probit',0,1);
%     sigmoidfit_TT('J190802c32','probit',0,1);
%     sigmoidfit_TT('J190722b42','probit',0,1);
%     sigmoidfit_TT('J190731c41','probit',0,1);

end

for imonkey = 1:nmonkeys
monkey_ana = monkeys_ana{imonkey};

try 
filename = ['pop_behav_ana_summary_',monkey_ana];
eval(['load ', filename])
% exnovo

catch
    
nsessions_all = 0;
nsessions_ana = 0;

%sessions
eval(['sessionlist_',monkey_ana]);

% steepness_all = [];
% rho_all = [];
% rhoABBA_all = [];
% hyst_all = [];
% %
% miu_ord_all = [];
% %
% steepness_SOABBA_sep_all = [];
% rho_SOABBA_sep_all = [];
% %
% quantrange_all = [];
% 
% RT_all = [];
% RT_JC_all = [];
% RT_SO_all = [];
% RT_AB_all = [];
% RT_BA_all = [];
% 
% RT_JCchA_all = [];
% RT_SOchA_all = [];
% RT_ABchA_all = [];
% RT_BAchA_all = [];
% RT_JCchB_all = [];
% RT_SOchB_all = [];
% RT_ABchB_all = [];
% RT_BAchB_all = [];
% 
% ERratio_all = [];
% percCor_all = [];
% fractofA_all = [];
% 
% sessionnames_ana_all = {};
% cellnames_ana_all = {};
% 
% juicetypes_all = {};
%    
% % ndynrange = [];
% nosatur = [];

meanOVA_all = []; % average value of offer A
meanOVB_all = []; % average value of offer B
meanCVA_all = []; % average value of chosen offer A
meanCVB_all = []; % average value of chosen offer B

for isession = 1:size(sessions,1)
    session = sessions{isession};
	readsession_TT

    ncells_isess = size(cells_td,1);
    ncells_tgtarea = 0;
    for icell = 1:ncells_isess  
        cellname = [num2str(cells_td(icell,1)),num2str(cells_td(icell,2))];
        if isequal(brainarea, arearead_TT(session,cellname)) %&& ~ismember([session,cellname],doubles)
            ncells_tgtarea = ncells_tgtarea + 1;
        end
    end
    
    if ncells_tgtarea > 0
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % remove sessions accordingly
       % skip sessions with too few trials
        filename = [dirroot,session,cellname,'_data'];
        eval(['load ',filename])
        clear celldata celldataerror trace trialRecord
        if size(goodTrials_JC,1)<mintrialnum | size(goodTrials_SO,1)<mintrialnum 
            continue
        end   
        nsessions_all = nsessions_all + 1;      
        
        if removesessions | ~removesessions
%             if psyphydata(find(psyphydata(:,3)==30,1),1) - psyphydata(find(psyphydata(:,3)==25,1),1) >700
%                 continue
%             end
            % remove sessions which has bad dynamic ranges -WS
            % bad dynamic ranges mean no saturation point at either end of sigmoidal curve
            filename = [dirroot,session,cellname,'_tuning'];
            eval(['load ',filename])
            table01_JC   = tuning.JC.AB.table01;
            table01_SOAB = tuning.SO.AB.table01;
            table01_SOBA = tuning.SO.BA.table01;
            table01_SO   = tuning.SO.ABA.table01;
            % re-organize based on quantity ratio
            types = {'JC','SOAB','SOBA','SO'};
            for itype = 1:4
                type = types{itype};
                eval(['table01 = table01_',type,';'])
                quantityratio = table01(:,1)./table01(:,2);
                qratios_uni = unique(quantityratio);
                nqratios = length(qratios_uni);
                ratiotable01 = zeros(nqratios,4);
                for iqratio = 1:nqratios
                    qratio_iqr = qratios_uni(iqratio);
                    ind_iqr = ismember(quantityratio,qratio_iqr);
                    table01_iqr = table01(ind_iqr,:);
                    ratiotable01(iqratio,1:2) = table01_iqr(1,1:2);
                    ratiotable01(iqratio,3)   = sum(table01_iqr(:,3).*table01_iqr(:,4))./sum(table01_iqr(:,4));
                    ratiotable01(iqratio,4)   = sum(table01_iqr(:,4)); 
                    % ratiotable01(iqratio,3)   = min(table01_iqr(:,3));
                    % ratiotable01(iqratio,4)   = min(table01_iqr(:,4)); 
                end
                eval(['ndynrange.',type,'(nsessions_all,:) = sum(ratiotable01(:,3)<0.9 & ratiotable01(:,3)>0.1);'])
            end
            %
%             nosatur.choiceB(nsessions_all,:) = min([max(table01_JC(1:end-1,3)),max(table01_SOAB(1:end-1,3)),max(table01_SOBA(1:end-1,3))])<0.9;
%             nosatur.choiceA(nsessions_all,:) = max([min(table01_JC(1:end-1,3)),min(table01_SOAB(1:end-1,3)),min(table01_SOBA(1:end-1,3))])>0.1;
            nosatur.choiceB(nsessions_all,:) = min([max(table01_JC(1:end-1,3)),max(table01_SO(1:end-1,3))])<0.9;
            nosatur.choiceA(nsessions_all,:) = max([min(table01_JC(1:end-1,3)),min(table01_SO(1:end-1,3))])>0.1;
            if removesessions & 0 %in the new version, remove these sessions later
            if ndynrange.JC(nsessions_all,:)<1 || ndynrange.SOAB(nsessions_all,:)<1 || ndynrange.SOBA(nsessions_all,:)<1
                disp(['session #',session])
                disp('improper dynamic range!')
                continue
            end
            if nosatur.choiceB(nsessions_all,:)
                disp(['session #',session])
                disp('improper dynamic range at choice B end!')
                continue
            end
            if nosatur.choiceA(nsessions_all,:)
                disp(['session #',session])
                disp('improper dynamic range at choice A end!')
                continue
            end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        nsessions_ana = nsessions_ana + 1; 
        sessionnames_ana_all{nsessions_ana,1} = session;
        cellnames_ana_all{nsessions_ana,1} = [session,cellname];
        
        
        % sigmoidal fitting
        % 
        warning off
        % [psyphycell] = sigmoidfit_TT_OrdChHyst([session,cellname],'probit','Only',1); % 'Only': hysteresis SO on SO or JC on JC; 'Both': SO and JC on either SO or JC
        % [psyphycell] = sigmoidfit_TT_OrdChHyst([session,cellname],'probit','Both',1); % 'Only': hysteresis SO on SO or JC on JC; 'Both': SO and JC on either SO or JC
%         [psyphycell] = sigmoidfit_TT_OrdChHyst([session,cellname],'logit','Both','log',1); % 'Only': hysteresis SO on SO or JC on JC; 'Both': SO and JC on either SO or JC
        [psyphycell] = sigmoidfit_TT_OrdChHyst([session,cellname],dologitorprobit(2:end),'Both','log',1); % 'Only': hysteresis SO on SO or JC on JC; 'Both': SO and JC on either SO or JC
        %
        % steepness
        if fittingwithCH
            steepness_JC = psyphycell.JC.OrdChHyst.sigmoidfit.beta(2);
            steepness_SO = psyphycell.SO.OrdChHyst.sigmoidfit.beta(2);
        elseif ~fittingwithCH
            steepness_JC = psyphycell.JC.NonChHyst.sigmoidfit.beta(2);
            steepness_SO = psyphycell.SO.NonChHyst.sigmoidfit.beta(2);  
        end
        %
        % relative value
        if fittingwithCH
            rho_JC = exp(-(psyphycell.JC.OrdChHyst.sigmoidfit.beta(1))/(psyphycell.JC.OrdChHyst.sigmoidfit.beta(2)));
            rho_SO = exp(-(psyphycell.SO.OrdChHyst.sigmoidfit.beta(1))/(psyphycell.SO.OrdChHyst.sigmoidfit.beta(2)));
            rho_SOAB = exp(-(psyphycell.SO.OrdChHyst.sigmoidfit.beta(1)-psyphycell.SO.OrdChHyst.sigmoidfit.beta(3))/(psyphycell.SO.OrdChHyst.sigmoidfit.beta(2)));
            rho_SOBA = exp(-(psyphycell.SO.OrdChHyst.sigmoidfit.beta(1)+psyphycell.SO.OrdChHyst.sigmoidfit.beta(3))/(psyphycell.SO.OrdChHyst.sigmoidfit.beta(2)));
        elseif ~fittingwithCH
            rho_JC = exp(-(psyphycell.JC.NonChHyst.sigmoidfit.beta(1))/(psyphycell.JC.NonChHyst.sigmoidfit.beta(2)));
            rho_SO = exp(-(psyphycell.SO.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SO.NonChHyst.sigmoidfit.beta(2)));
            rho_SOAB = exp(-(psyphycell.SO.NonChHyst.sigmoidfit.beta(1)-psyphycell.SO.NonChHyst.sigmoidfit.beta(3))/(psyphycell.SO.NonChHyst.sigmoidfit.beta(2)));
            rho_SOBA = exp(-(psyphycell.SO.NonChHyst.sigmoidfit.beta(1)+psyphycell.SO.NonChHyst.sigmoidfit.beta(3))/(psyphycell.SO.NonChHyst.sigmoidfit.beta(2)));
        end
        %
        % one step choice hysteresis (both on SO and both on JC)
        hyst_JC = (psyphycell.JC.OrdChHyst.sigmoidfit.beta(3))/(psyphycell.JC.OrdChHyst.sigmoidfit.beta(2));
        hyst_SO = (psyphycell.SO.OrdChHyst.sigmoidfit.beta(4))/(psyphycell.SO.OrdChHyst.sigmoidfit.beta(2));
        %
        % order bias
        if fittingwithCH
            miu_ord_SO = -(psyphycell.SO.OrdChHyst.sigmoidfit.beta(3))/(psyphycell.SO.OrdChHyst.sigmoidfit.beta(2));
        elseif ~fittingwithCH
            miu_ord_SO = -(psyphycell.SO.NonChHyst.sigmoidfit.beta(3))/(psyphycell.SO.NonChHyst.sigmoidfit.beta(2));
        end
        %
        % SO steepness and rho, separate for AB and BA
        if fittingwithCH
%             steepness_SOAB_sep = psyphycell.SOAB.OrdChHyst.sigmoidfit.beta(2);
%             steepness_SOBA_sep = psyphycell.SOBA.OrdChHyst.sigmoidfit.beta(2);
%             rho_SOAB_sep = exp(-(psyphycell.SOAB.OrdChHyst.sigmoidfit.beta(1))/(psyphycell.SOAB.OrdChHyst.sigmoidfit.beta(2)));
%             rho_SOBA_sep = exp(-(psyphycell.SOBA.OrdChHyst.sigmoidfit.beta(1))/(psyphycell.SOBA.OrdChHyst.sigmoidfit.beta(2)));
            steepness_SOAB_sep = psyphycell.SOAB.NonChHyst.sigmoidfit.beta(2);
            steepness_SOBA_sep = psyphycell.SOBA.NonChHyst.sigmoidfit.beta(2);
            rho_SOAB_sep = exp(-(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(2)));
            rho_SOBA_sep = exp(-(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(2)));
        elseif ~fittingwithCH
            steepness_SOAB_sep = psyphycell.SOAB.NonChHyst.sigmoidfit.beta(2);
            steepness_SOBA_sep = psyphycell.SOBA.NonChHyst.sigmoidfit.beta(2);
            rho_SOAB_sep = exp(-(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(2)));
            rho_SOBA_sep = exp(-(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(2)));
        end
        
        % %
        % % 
        % flip if needed
        if (rho_JC + rho_SO)/2 < 1
        % if (rho_JC + rho_SOAB + rho_SOBA)/3 < 1
        % if rho_JC < 1
        % if min([rho_JC, rho_SOAB, rho_SOBA]) < 1
            rho_JC = 1/rho_JC;
            rho_SO = 1/rho_SO;
            rho_SOAB_sep = 1/rho_SOAB_sep;
            rho_SOBA_sep = 1/rho_SOBA_sep;
            rho_SOAB = 1/rho_SOAB;
            rho_SOBA = 1/rho_SOBA;
        end
        
        steepness_all(nsessions_ana,:) = [steepness_JC,steepness_SO];
        rho_all(nsessions_ana,:) = [rho_JC,rho_SO];
        rhoABBA_all(nsessions_ana,:) = [rho_SOAB,rho_SOBA];           
        miu_ord_all(nsessions_ana,:) = miu_ord_SO;
        hyst_all(nsessions_ana,:) = [hyst_JC, hyst_SO];
        steepness_SOABBA_sep_all(nsessions_ana,:) = [steepness_SOAB_sep, steepness_SOBA_sep];
        rho_SOABBA_sep_all(nsessions_ana,:) = [rho_SOAB_sep, rho_SOBA_sep];    
        
        % % %
        % fraction of choosing A
        filename = [dirroot,session,cellname,'_tuning'];
        eval(['load ',filename])
        neuract_JC = tuning.JC.AB.neuract.bytrial.preoffer;
        chosenJ_JC = neuract_JC(:,4).*neuract_JC(:,5);
        fractofA_JC = sum(chosenJ_JC==1)./length(chosenJ_JC);
        neuract_SO = tuning.SO.ABA.neuract.bytrial.preoffers;
        chosenJ_SO = neuract_SO(:,5);
        fractofA_SO = sum(chosenJ_SO==1)./length(chosenJ_SO);
        fractofA_all(nsessions_ana,:) = [fractofA_JC,fractofA_SO];
        
        % % %
        % meanOVA meanOVB meanCVA meanCVB        
        if (rho_JC + rho_SO)/2 < 1
            OVA_JC = neuract_JC(:,2)./rho_JC;   
            OVA_SO = neuract_SO(:,2)./rho_SO;   
        else
            OVA_JC = neuract_JC(:,2).*rho_JC;
            OVA_SO = neuract_SO(:,2).*rho_SO;
        end
        OVB_JC = neuract_JC(:,3);
        OVB_SO = neuract_SO(:,3);
        CVA_JC = OVA_JC;
        CVA_JC(chosenJ_JC == -1) = [];
        CVB_JC = OVB_JC;
        CVB_JC(chosenJ_JC ==  1) = [];
        CVA_SO = OVA_SO;
        CVA_SO(chosenJ_SO == -1) = [];
        CVB_SO = OVB_SO;
        CVB_SO(chosenJ_SO ==  1) = [];
        meanOVA_all(nsessions_ana,:) = [nanmean(OVA_JC), nanmean(OVA_SO)];
        meanOVB_all(nsessions_ana,:) = [nanmean(OVB_JC), nanmean(OVB_SO)];
        meanCVA_all(nsessions_ana,:) = [nanmean(CVA_JC), nanmean(CVA_SO)];
        meanCVB_all(nsessions_ana,:) = [nanmean(CVB_JC), nanmean(CVB_SO)];
        
        
        % % %
        % quantity range and juice type
        filename = [dirroot,session,cellname,'_data'];
        eval(['load ',filename])
        clear celldata celldataerror trace
        %
        quantrangeA = [min(abs(goodTrials_JC(:,2))), max(abs(goodTrials_JC(:,2)))];
        quantrangeB = [min(abs(goodTrials_JC(:,3))), max(abs(goodTrials_JC(:,3)))];
        valuerangeA = quantrangeA;
        valuerangeB = quantrangeB;
        quantrange_all(nsessions_ana,:) = [valuerangeA(2)-valuerangeA(1),valuerangeB(2)-valuerangeB(1)];
        %
        % juice type
        for ind_1stnonforce = 1:size(trialRecord,2)
            ind_goodorder = [trialRecord(ind_1stnonforce).currentOffers.goodId];
            ind_findzero = find(ind_goodorder==0);
            if isempty(ind_findzero)
                break
            end
        end        
        ind_goodorder = [trialRecord(ind_1stnonforce).currentOffers.goodId];
        goodA = trialRecord(ind_1stnonforce).currentOffers(ind_goodorder(1)).good.juice;
        goodB = trialRecord(ind_1stnonforce).currentOffers(ind_goodorder(2)).good.juice;
        juicetypes_all(nsessions_ana,:) = {goodA, goodB};
        
        
        % Reaction time and error rate 
        %
        % reaction time
        if doreacttime
        filename = [dirroot,session,cellname,'_data'];
        eval(['load ',filename])
        clear celldata celldataerror trace trialRecord
        %
        RT_JC = [];
        trial_JC = goodTrials_JC(:,1);
        for itrial = 1:length(trial_JC)
            num_itrial = trial_JC(itrial);
            psyphydata_itrial = psyphydata(psyphydata(:,2)==num_itrial,:);
            try
                RT_JC(itrial) = psyphydata_itrial(psyphydata_itrial(:,3)==39,1)-psyphydata_itrial(psyphydata_itrial(:,3)==35,1);
                if RT_JC(itrial)<=0
                    RT_JC(itrial) = nan;
                end
            catch
                RT_JC(itrial) = nan;
            end
        end
        %
%         RT_SO = [];
        trial_SO = goodTrials_SO(:,1);
        for itrial = 1:length(trial_SO)
            num_itrial = trial_SO(itrial);
            psyphydata_itrial = psyphydata(psyphydata(:,2)==num_itrial,:);
            try
                RT_SO(itrial) = psyphydata_itrial(psyphydata_itrial(:,3)==39,1)-psyphydata_itrial(psyphydata_itrial(:,3)==27,1);
                if RT_SO(itrial)<=0
                    RT_SO(itrial) = nan;
                end
            catch
                RT_SO(itrial) = nan;
            end
        end
        filename = [dirroot,session,cellname,'_tuning'];
        eval(['load ', filename])
        ABBAs = tuning.SO.ABA.neuract.bytrial.preoffers(:,[1,6]);
        [~,iii] = ismember(goodTrials_SO(:,1),ABBAs(:,1));
        ind_ABBAs = ABBAs(iii,2);
        RT_JC_all(nsessions_ana,:) = nanmean(RT_JC);
        RT_SO_all(nsessions_ana,:) = nanmean(RT_SO);
        RT_AB_all(nsessions_ana,:) = nanmean(RT_SO(ind_ABBAs== 1));
        RT_BA_all(nsessions_ana,:) = nanmean(RT_SO(ind_ABBAs==-1)); 
        %
        chA_JC = [tuning.JC.AB.neuract.bytrial.preoffer(:,1), ...
                  tuning.JC.AB.neuract.bytrial.preoffer(:,4).*tuning.JC.AB.neuract.bytrial.preoffer(:,5)]; % 1:A; -1:B
        [~,iii] = ismember(goodTrials_JC(:,1),chA_JC(:,1));
        ind_chA_JC = chA_JC(iii,2);      
        RT_JCchA_all(nsessions_ana,:) = nanmean(RT_JC(ind_chA_JC== 1));
        RT_JCchB_all(nsessions_ana,:) = nanmean(RT_JC(ind_chA_JC==-1));
        % 
        chA_SO = [tuning.SO.ABA.neuract.bytrial.preoffers(:,1), tuning.SO.ABA.neuract.bytrial.preoffers(:,5)]; % 1:A; -1:B
        [~,iii] = ismember(goodTrials_SO(:,1),chA_SO(:,1));
        ind_chA_SO = chA_SO(iii,2); 
        RT_SOchA_all(nsessions_ana,:) = nanmean(RT_SO(ind_chA_SO== 1));
        RT_SOchB_all(nsessions_ana,:) = nanmean(RT_SO(ind_chA_SO==-1));
        RT_ABchA_all(nsessions_ana,:) = nanmean(RT_SO(ind_ABBAs== 1 & ind_chA_SO== 1));
        RT_ABchB_all(nsessions_ana,:) = nanmean(RT_SO(ind_ABBAs== 1 & ind_chA_SO==-1));
        RT_BAchA_all(nsessions_ana,:) = nanmean(RT_SO(ind_ABBAs==-1 & ind_chA_SO== 1));
        RT_BAchB_all(nsessions_ana,:) = nanmean(RT_SO(ind_ABBAs==-1 & ind_chA_SO==-1));  
        else
        RT_JC_all = [];
        RT_SO_all = [];
        RT_AB_all = [];
        RT_BA_all = [];
        RT_JCchA_all = [];
        RT_JCchB_all = [];
        RT_SOchA_all = [];
        RT_SOchB_all = [];
        RT_ABchA_all = [];
        RT_ABchB_all = [];
        RT_BAchA_all = [];
        RT_BAchB_all = [];
        end        
        %
                          
    end % isequal brainarea
end% for isession

% save files
filesavename = ['pop_behav_ana_summary_',monkey_ana];
eval(['save ', filesavename,' ndynrange nosatur nsessions_all nsessions_ana sessionnames_ana_all cellnames_ana_all juicetypes_all '...
      'steepness_all rho_all rhoABBA_all hyst_all miu_ord_all rho_SOABBA_sep_all steepness_SOABBA_sep_all quantrange_all '...
      'RT_JC_all RT_SO_all RT_AB_all RT_BA_all ERratio_all percCor_all fractofA_all '...
      'RT_JCchA_all RT_SOchA_all RT_ABchA_all RT_BAchA_all RT_JCchB_all RT_SOchB_all RT_ABchB_all RT_BAchB_all '...
      'meanOVA_all meanOVB_all meanCVA_all meanCVB_all ' ])

end


examplelist = ismember(sessionnames_ana_all,examplesessions);

% remove sessions with incorrect steepness
if remove_steepnessout & ~removesessions
    if isequal(monkey_ana,'both')
        monkeynames_all = [];
        for isession = 1:nsessions_all
            monkeynames_all(isession,1) = sessionnames_ana_all{isession}(1);           
        end
        %
        % Juan
        ind_J = ismember(monkeynames_all,'J');
        %
        quantiles_JC = quantile(steepness_all(ind_J,1),[0.25 0.5 0.75]);
        IQR_JC = quantiles_JC(3) - quantiles_JC(1);
        steepoutlier_JC = [quantiles_JC(1)-1.5*IQR_JC, quantiles_JC(3)+1.5*IQR_JC];
        %
        quantiles_SO = quantile(steepness_all(ind_J,2),[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-1.5*IQR_SO, quantiles_SO(3)+1.5*IQR_SO];
        %
        deltasteep = steepness_all(ind_J,1) - steepness_all(ind_J,2);
        quantiles_del = quantile(deltasteep,[0.25 0.5 0.75]);
        IQR_del = quantiles_del(3) - quantiles_del(1);
        steepoutlier_del = [quantiles_del(1)-1.5*IQR_del, quantiles_del(3)+1.5*IQR_del];
        % 
        ind_goodJC = steepness_all(ind_J,1) > steepoutlier_JC(1) & steepness_all(ind_J,1) < steepoutlier_JC(2);
        ind_goodSO = steepness_all(ind_J,2) > steepoutlier_SO(1) & steepness_all(ind_J,2) < steepoutlier_SO(2);
        ind_goodAB = steepness_SOABBA_sep_all(ind_J,1) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(ind_J,1) < steepoutlier_SO(2);
        ind_goodBA = steepness_SOABBA_sep_all(ind_J,2) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(ind_J,2) < steepoutlier_SO(2);
        % ind_good = ind_goodJC & ind_goodAB & ind_goodBA;
        ind_good_J = ind_goodJC & ind_goodSO;
        % ind_good = ind_goodJC & ind_goodSO & ind_goodAB & ind_goodBA;
        nsession_steepnessout_J = size(ind_good_J,1)-sum(ind_good_J);
        %
        % Gervinho
        ind_G = ismember(monkeynames_all,'G');
        %
        quantiles_JC = quantile(steepness_all(ind_G,1),[0.25 0.5 0.75]);
        IQR_JC = quantiles_JC(3) - quantiles_JC(1);
        steepoutlier_JC = [quantiles_JC(1)-1.5*IQR_JC, quantiles_JC(3)+1.5*IQR_JC];
        %
        quantiles_SO = quantile(steepness_all(ind_G,2),[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-1.5*IQR_SO, quantiles_SO(3)+1.5*IQR_SO];
        %
        deltasteep = steepness_all(ind_G,1) - steepness_all(ind_G,2);
        quantiles_del = quantile(deltasteep,[0.25 0.5 0.75]);
        IQR_del = quantiles_del(3) - quantiles_del(1);
        steepoutlier_del = [quantiles_del(1)-1.5*IQR_del, quantiles_del(3)+1.5*IQR_del];
        % 
        ind_goodJC = steepness_all(ind_G,1) > steepoutlier_JC(1) & steepness_all(ind_G,1) < steepoutlier_JC(2);
        ind_goodSO = steepness_all(ind_G,2) > steepoutlier_SO(1) & steepness_all(ind_G,2) < steepoutlier_SO(2);
        ind_goodAB = steepness_SOABBA_sep_all(ind_G,1) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(ind_G,1) < steepoutlier_SO(2);
        ind_goodBA = steepness_SOABBA_sep_all(ind_G,2) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(ind_G,2) < steepoutlier_SO(2);
        % ind_good = ind_goodJC & ind_goodAB & ind_goodBA;
        ind_good_G = ind_goodJC & ind_goodSO;
        % ind_good = ind_goodJC & ind_goodSO & ind_goodAB & ind_goodBA;
        nsession_steepnessout_G = size(ind_good_G,1)-sum(ind_good_G);
        %       
        ind_good = [ind_good_G; ind_good_J];
        nsession_steepnessout = nsession_steepnessout_G + nsession_steepnessout_J;
    else
        quantiles_JC = quantile(steepness_all(:,1),[0.25 0.5 0.75]);
        IQR_JC = quantiles_JC(3) - quantiles_JC(1);
        steepoutlier_JC = [quantiles_JC(1)-1.5*IQR_JC, quantiles_JC(3)+1.5*IQR_JC];
        %
        quantiles_SO = quantile(steepness_all(:,2),[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-1.5*IQR_SO, quantiles_SO(3)+1.5*IQR_SO];
        %
        deltasteep = steepness_all(:,1) - steepness_all(:,2);
        quantiles_del = quantile(deltasteep,[0.25 0.5 0.75]);
        IQR_del = quantiles_del(3) - quantiles_del(1);
        steepoutlier_del = [quantiles_del(1)-1.5*IQR_del, quantiles_del(3)+1.5*IQR_del];
        % 
        ind_goodJC = steepness_all(:,1) > steepoutlier_JC(1) & steepness_all(:,1) < steepoutlier_JC(2);
        ind_goodSO = steepness_all(:,2) > steepoutlier_SO(1) & steepness_all(:,2) < steepoutlier_SO(2);
        ind_goodAB = steepness_SOABBA_sep_all(:,1) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(:,1) < steepoutlier_SO(2);
        ind_goodBA = steepness_SOABBA_sep_all(:,2) > steepoutlier_SO(1) & steepness_SOABBA_sep_all(:,2) < steepoutlier_SO(2);
    %     ind_good = ind_goodJC & ind_goodAB & ind_goodBA;
        ind_good = ind_goodJC & ind_goodSO;
        % ind_good = ind_goodJC & ind_goodSO & ind_goodAB & ind_goodBA;
        nsession_steepnessout = size(ind_good,1)-sum(ind_good);    
    end
    
elseif ~remove_steepnessout & ~removesessions
    nsession_steepnessout = 0;
    ind_good = logical(ones(length(steepness_all),1));
end



% removesessions: based on dynamic range and saturation
%                 used to be done earlier, in the new version, remove these
%                 sessions here
if removesessions & ~remove_steepnessout
    ind_good = ndynrange.JC > 1 & ndynrange.SO > 1;
    nsession_steepnessout = size(ind_good,1)-sum(ind_good);
elseif ~remove_steepnessout & ~removesessions
    nsession_steepnessout = 0;
    ind_good = logical(ones(length(steepness_all),1));
end


% % % 
% plot
% % % 
% steepness and rho and order bias
% basic plots
if onlyplotpapers | ~onlyplotpapers
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot = length(htt);
    else
        figure(isumplot);
    end    
    % sigmoidal curve
    subplot(nmonkeys,4,1+(imonkey-1)*4);
    cellname = examplecells{imonkey};
    plot_sigmoidcurve(cellname);
    % rho
    subplot(nmonkeys,4,2+(imonkey-1)*4);
    % plot(rho_all(ind_good,1),rho_all(ind_good,2),'ko','MarkerSize',9,'LineWidth',1);hold on
    XX = [0 7.5]; YY = [0 7.5];
    % XX = [0 ceil(max(max(rho_all(ind_good,:))))]; YY = XX;
    hold on; plot(XX, YY, '--','Color', [0 0 0],'LineWidth',1);
    Sigma_ell = cov(rho_all(ind_good,1),rho_all(ind_good,2));
    mu_ell(1) = mean(rho_all(ind_good,1));
    mu_ell(2) = mean(rho_all(ind_good,2));
    hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90)
    plot(rho_all(ind_good,1),rho_all(ind_good,2),'ko','MarkerSize',9,'LineWidth',1);hold on
    try
    sscat1 = scatter(rho_all(logical(examplelist),1),rho_all(logical(examplelist),2),...
                     90,'o','markerfacecolor',[1 0.5 0],'markeredgecolor','k');
    sscat1.MarkerFaceAlpha = .6;
    end
    [~, pp_ttest] = ttest(rho_all(ind_good,1),rho_all(ind_good,2));
    [pp_wil, ~  ] = signrank(rho_all(ind_good,1),rho_all(ind_good,2));
    [rrr,p_rrr] = corr(rho_all(ind_good,1),rho_all(ind_good,2));
%     text(XX(1)+XX(2)/50,YY(2)-YY(2)/8,{['corr = ',num2str(rrr,'%0.2f')];...
%                                        ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
%                                        ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
%                                        ['mean \Delta = ',num2str(mean(-rho_all(ind_good,1)+rho_all(ind_good,2)),'%.2f')];...
%                                        ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']},'fontsize',11);
    text(XX(1)+XX(2)/50,YY(2)-YY(2)/8,{['mean \Delta = ',num2str(mean(-rho_all(ind_good,1)+rho_all(ind_good,2)),'%.2f')];...
                                       ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                       ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                      },'fontsize',11);
    xlabel('\rho Task 1'); ylabel('\rho Task 2');
    % title({['monkey ',monkey_ana]});
    axis([XX YY]); axis square; box off
    xTickrange = round(XX(1)):round((round(XX(2))-round(XX(1)))/3):round(XX(2));
    xTickLabelrange = num2str(xTickrange');
    yTickrange = round(YY(1)):round((round(YY(2))-round(YY(1)))/3):round(YY(2));
    yTickLabelrange = num2str(yTickrange');
    set(gca, 'XTick', xTickrange, ...                  
             'XTickLabel', xTickLabelrange, ...                    
             'YTick', yTickrange, ...
             'YTickLabel', yTickLabelrange);
    set(gca,'fontsize',13)    
    % steepness
    subplot(nmonkeys,4,3+(imonkey-1)*4);
    % plot(steepness_all(ind_good,1),steepness_all(ind_good,2),'ko','MarkerSize',9,'LineWidth',1); hold on
    XX = [0 13]; YY = [0 13];
    % XX = [0 ceil(max(max(steepness_all(ind_good,:))))]; YY = XX;
    hold on; plot(XX, YY, '--','Color', [0 0 0],'LineWidth',1);
    Sigma_ell = cov(steepness_all(ind_good,1),steepness_all(ind_good,2));
    mu_ell(1) = mean(steepness_all(ind_good,1));
    mu_ell(2) = mean(steepness_all(ind_good,2));
    hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
    plot(steepness_all(ind_good,1),steepness_all(ind_good,2),'ko','MarkerSize',9,'LineWidth',1); hold on
    try
    % plot(steepness_all(logical(examplelist),1),steepness_all(logical(examplelist),2),'ko','MarkerFaceColor','b','MarkerSize',9,'LineWidth',1);
    sscat1 = scatter(steepness_all(logical(examplelist),1),steepness_all(logical(examplelist),2),...
                     90,'o','markerfacecolor',[1 0.5 0],'markeredgecolor','k');
    sscat1.MarkerFaceAlpha = .6;
    end
    [~, pp_ttest] = ttest(steepness_all(ind_good,1),steepness_all(ind_good,2));
    [pp_wil, ~  ] = signrank(steepness_all(ind_good,1),steepness_all(ind_good,2));
    text(XX(1)+XX(2)/50,YY(2)-YY(2)/8,{['mean \Delta = ',num2str(mean(-steepness_all(ind_good,1)+steepness_all(ind_good,2)),'%.2f')];...
                                       ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                       ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                      },'fontsize',11);
    xlabel('\eta Task 1'); ylabel('\eta Task 2');
    % title({['monkey ',monkey_ana]});
    axis([XX YY]); axis square; box off
    xTickrange = round(XX(1)):round((round(XX(2))-round(XX(1)))/3):round(XX(2));
    xTickLabelrange = num2str(xTickrange');
    yTickrange = round(YY(1)):round((round(YY(2))-round(YY(1)))/3):round(YY(2));
    yTickLabelrange = num2str(yTickrange');
    set(gca, 'XTick', xTickrange, ...                  
             'XTickLabel', xTickLabelrange, ...                    
             'YTick', yTickrange, ...
             'YTickLabel', yTickLabelrange);
    set(gca,'fontsize',13)         
    % miu_ord
    subplot(nmonkeys,4,4+(imonkey-1)*4);
    XX = [-1, 2]; 
    if isequal(monkey_ana,'Gervinho')
        YY = [0, 38];
    elseif isequal(monkey_ana,'Juan')
        YY = [0, 21];
    elseif isequal(monkey_ana,'both')
        YY = [0, 60];
    end
    edges = [XX(1):0.1:XX(2)];
    % histogram(miu_ord_all(ind_good),edges,'FaceColor',[.4 .4 .4]); hold on;
    % plot([0 0],YY, '--','Color', [0 0 0],'LineWidth',1); hold on;
    histogram(2*miu_ord_all(ind_good).*rho_all(ind_good,2),edges,'FaceColor',[.4 .4 .4]); hold on;
    % plot(2*miu_ord_all(logical(examplelist)).*rho_all(logical(examplelist),2), YY(2)*0.75, 'v', 'MarkerSize', 6,'MarkerFaceColor','b','MarkerEdgeColor','k');
    plot(mean(2*miu_ord_all(ind_good).*rho_all(ind_good,2)), YY(2)*0.75, '^', 'MarkerSize', 6,'MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor',[1 0.5 0]);
    [~, pp_ttest] = ttest(2*miu_ord_all(ind_good).*rho_all(ind_good,2));
    [ pp_wil,  ~] = signrank(2*miu_ord_all(ind_good).*rho_all(ind_good,2));
    text(XX(1)+XX(2)/1.5,YY(2)-YY(2)/8, {['mean = ',num2str(mean(2*miu_ord_all(ind_good).*rho_all(ind_good,2)),'%.2f')];...
                                         ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                         ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                        },'fontsize',11);
    xlabel('\epsilon Task 2');
    ylabel('N sessions');
    % title({['monkey ',monkey_ana]});
    axis([XX YY]); axis square; box off
    set(gca,'TickDir','out');
    set(gca,'fontsize',13);   
    %
    % add monkey names
    if isequal(monkey_ana,monkeys_ana{1})
        axes('position',[.055 .70 .2 .05]);
    elseif isequal(monkey_ana,monkeys_ana{2})
        axes('position',[.055 .25 .2 .05]); 
    end
    hold on; h = text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    set(h,'Rotation',90);
end


% % %
% rho JC AB and BA
% rho AB BA JC
if onlyplotpapers
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1650 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot35 = length(htt);
    else
        figure(isumplot35);
    end
    rho_JC_all = rho_all(:,1);
    rho_SO_all = rho_all(:,2);
    rho_AB_all = rhoABBA_all(:,1);
    rho_BA_all = rhoABBA_all(:,2);
    plottypes1 = {'JC','JC','JC'};
    plottypes2 = {'SO','AB','BA'};
%     xlabelnames = {'log relative value \rho of Task 1', 'log relative value \rho of Task 1',     'log relative value \rho of Task 1',   };
%     ylabelnames = {'log relative value \rho of Task 2', 'log relative value \rho of Task 2(AB)', 'log relative value \rho of Task 2(BA)'};
    xlabelnames = {'relative value \rho of Task 1', 'relative value \rho of Task 1',     'relative value \rho of Task 1',   };
    ylabelnames = {'relative value \rho of Task 2', 'relative value \rho of Task 2(AB)', 'relative value \rho of Task 2(BA)'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes+1,iplottype+(imonkey-1)*(nplottypes+1));
        eval(['XXX = rho_',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = rho_',plottypename2,'_all(ind_good,:);'])
        %
        XXX_lim = [rho_JC_all(ind_good,:); rho_AB_all(ind_good,:); rho_BA_all(ind_good,:)];       
        YYY_lim = XXX_lim;   
        %
        % plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        aa = deming(XXX,YYY);
        % XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5];
        % YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5];
%         XX = [0,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20];
%         YY = [0,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/20];
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/20,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/20];
        YYfit = aa(2)*XX+aa(1);
        [~, pp_ttest] = ttest(XXX, YYY);
        [pp_wil, ~  ] = signrank(XXX, YYY);
        hold on; plot(XX,YY,'--','LineWidth',1,'color',[0.7 0.7 0.7]);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        plot(XXX, YYY, 'ko','MarkerSize',6,'LineWidth',1);
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/8,...
            {['t test: p = ',num2str(pp_ttest,'%1.1g')]; ...
             ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')]; ...             
             ['mean \Delta = ',num2str(mean(-XXX+YYY),'%.2f')]; ...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions'];...
             }, 'fontsize', 9);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
%         title({[monkey_ana];['analyzed session # =',num2str(sum(ind_good))]});
        % xTickrange = round(XX(1),0):round((round(XX(2),0)-round(XX(1),0))/3):round(XX(2),0);
%         xTickrange = round(XX(1),0):2:round(XX(2),0);
%         xTickLabelrange = num2str(xTickrange');
%         % yTickrange = round(YY(1),0):round((round(YY(2),0)-round(YY(1),0))/3):round(YY(2),0);
%         yTickrange = round(YY(1),0):2:round(YY(2),0);
%         yTickLabelrange = num2str(yTickrange');
%         set(gca, 'XTick', xTickrange, ...                  
%                  'XTickLabel', xTickLabelrange, ...                    
%                  'YTick', yTickrange, ...
%                  'YTickLabel', yTickLabelrange);
        set(gca,'fontsize',11)    
        box off; axis([XX,YY]); axis square  
    end    
    %
    % plot rho v.s. rho difference
    rho_JC_all = (rho_all(:,1));
    rho_SO_all = (rho_all(:,2));
    rho_AB_all = (rhoABBA_all(:,1));
    rho_BA_all = (rhoABBA_all(:,2));
    rho_delta_all = (rho_all(:,2)) - (rho_all(:,1));
    rho_deltaAB_all = (rhoABBA_all(:,1)) - (rho_all(:,1));
    rho_deltaBA_all = (rhoABBA_all(:,2)) - (rho_all(:,1));
    steep_delta_all = steepness_all(:,2) - steepness_all(:,1);
    plottypes1 = {'rho_JC'};
    plottypes2 = {'rho_delta'};
    xlabelnames = {'relative value \rho of Task 1'};
    ylabelnames = {'\Delta relative value \rho (Task 2 - Task 1)'};
    nplottypes2 = length(plottypes1);
    nplottypes2 = 0;
    iplottype = 0;
    % % % 
    % rotation fitting (two formula)
    YYY_AB  = [rho_AB_all(ind_good,:)]-1; % rho_AB-1
    YYY_BA  = [rho_BA_all(ind_good,:)]-1; % rho_BA-1
    XXX1_AB = [ones(size(rho_AB_all(ind_good,:)))];
    XXX1_BA = [ones(size(rho_BA_all(ind_good,:)))];
    XXX2_JC = [rho_JC_all(ind_good,:)]-1; % rho_JC-1
    %
    mdl_AB = fitlm([XXX1_AB, XXX2_JC], YYY_AB,'Intercept',false);
    mdl_BA = fitlm([XXX1_BA, XXX2_JC], YYY_BA,'Intercept',false);
    aAB = mdl_AB.Coefficients.Estimate(1);
    aBA = mdl_BA.Coefficients.Estimate(1);
    bAB = mdl_AB.Coefficients.Estimate(2);
    bBA = mdl_BA.Coefficients.Estimate(2);
    paAB = coefTest(mdl_AB,[1,0],0);
    paBA = coefTest(mdl_BA,[1,0],0);
    pbAB = coefTest(mdl_AB,[0,1],1);
    pbBA = coefTest(mdl_BA,[0,1],1);
    %
    CI95_AB = coefCI(mdl_AB,0.05);
    CIaAB = CI95_AB(1,:);
    CIaABmean = mean(CIaAB);
    CIaABstd = CIaAB(2)-CIaABmean;
    CIbAB = CI95_AB(2,:);
    CIbABmean = mean(CIbAB);
    CIbABstd = CIbAB(2)-CIbABmean;
    %
    CI95_BA = coefCI(mdl_BA,0.05);
    CIaBA = CI95_BA(1,:);
    CIaBAmean = mean(CIaBA);
    CIaBAstd = CIaBA(2)-CIaBAmean;
    CIbBA = CI95_BA(2,:);
    CIbBAmean = mean(CIbBA);
    CIbBAstd = CIbBA(2)-CIbBAmean;
    %
    subplot(nmonkeys,nplottypes+nplottypes2+1,iplottype+nplottypes+1+(imonkey-1)*(nplottypes+nplottypes2+1));
    hold on; 
    h = text(0,0.5,{['\rhoTask2 - 1 = c1 + c2 * (\rhoTask1 - 1)'];...
                    ['AB: c1 = ',num2str(CIaABmean,'%.2f'),'\pm',num2str(CIaABstd,'%.2f'),'; p = ',num2str(paAB,'%1.1g')];...
                    ['AB: c2 = ',num2str(CIbABmean,'%.2f'),'\pm',num2str(CIbABstd,'%.2f'),'; p = ',num2str(pbAB,'%1.1g')];...
                    ['BA: c1 = ',num2str(CIaBAmean,'%.2f'),'\pm',num2str(CIaBAstd,'%.2f'),'; p = ',num2str(paBA,'%1.1g')];...
                    ['BA: c2 = ',num2str(CIbBAmean,'%.2f'),'\pm',num2str(CIbBAstd,'%.2f'),'; p = ',num2str(pbBA,'%1.1g')];...
                    },'FontSize',10);
    axis off
    %
    % add monkey names
    if isequal(monkey_ana,monkeys_ana{1})
        axes('position',[.055 .70 .2 .05]);
    elseif isequal(monkey_ana,monkeys_ana{2})
        axes('position',[.055 .25 .2 .05]); 
    end
    hold on; h = text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    set(h,'Rotation',90);
end


% % % 
% % %
% rho JC AB and BA
% rho AB BA vs JC ancova
if onlyplotpapers
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1650 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot36 = length(htt);
    else
        figure(isumplot36);
    end
    rho_JC_all = rho_all(:,1);
    rho_SO_all = rho_all(:,2);
    rho_AB_all = rhoABBA_all(:,1);
    rho_BA_all = rhoABBA_all(:,2);
    plottypes1 = {'JC','JC'};
    plottypes2 = {'AB','BA'};
%     xlabelnames = {'log relative value \rho of Task 1', 'log relative value \rho of Task 1',     'log relative value \rho of Task 1',   };
%     ylabelnames = {'log relative value \rho of Task 2', 'log relative value \rho of Task 2(AB)', 'log relative value \rho of Task 2(BA)'};
    xlabelnames = {'relative value \rho of Task 1',     'relative value \rho of Task 1',   };
    ylabelnames = {'relative value \rho of Task 2(AB)', 'relative value \rho of Task 2(BA)'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes+1,iplottype+(imonkey-1)*(nplottypes+1));
        eval(['XXX = rho_',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = rho_',plottypename2,'_all(ind_good,:);'])
        %
        XXX_lim = [rho_JC_all(ind_good,:); rho_AB_all(ind_good,:); rho_BA_all(ind_good,:)];       
        YYY_lim = XXX_lim;   
        aa = deming(XXX,YYY);
        % XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20];
        XX = [0 7.5];
        % YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/20,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/20];
        YY = [0 7.5];
        YYfit = aa(2)*XX+aa(1);
        [~, pp_ttest] = ttest(XXX, YYY);
        [pp_wil, ~  ] = signrank(XXX, YYY);
        hold on; plot(XX,YY,'--','LineWidth',1,'color',[0.7 0.7 0.7]);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        plot(XXX, YYY, 'ko','MarkerSize',6,'LineWidth',1);
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/8,...
            {['t test: p = ',num2str(pp_ttest,'%1.1g')]; ...
             ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')]; ...             
             ['mean \Delta = ',num2str(mean(-XXX+YYY),'%.2f')]; ...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions'];...
             }, 'fontsize', 9);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        set(gca,'fontsize',11)    
        box off; axis([XX,YY]); axis square  
    end
    %
    %
    % rotation fitting (ancova)
    XXX = [rho_JC_all(ind_good,:); rho_JC_all(ind_good,:)]; % rho_JC
    YYY = [rho_AB_all(ind_good,:); rho_BA_all(ind_good,:)]; % rho_SO
    ABBAgroup = [repelem({'AB'},length((rho_AB_all(ind_good,:)))),repelem({'BA'},length((rho_BA_all(ind_good,:))))]';
    % 
    % ANCOVA
    [h,atab,ctab,stats] = aoctool(XXX,YYY,ABBAgroup,0.05,'','','','off','parallel lines');
    % [h,atab,ctab,stats] = aoctool(XXX,YYY,ABBAgroup,0.05,'','','','off','separate lines');
    p_intc = ctab{3,5};
    try
        p_slpe = ctab{5,5};
    end
    %
    % fit two parellel lines
    YYY_AB  = [rho_AB_all(ind_good,:)]; % rho_AB
    YYY_BA  = [rho_BA_all(ind_good,:)]; % rho_BA
    YYY_ABBA = [YYY_AB; YYY_BA];
    XXX1_AB = [ones(size(rho_AB_all(ind_good,:)));  zeros(size(rho_BA_all(ind_good,:)))];
    XXX1_BA = [zeros(size(rho_AB_all(ind_good,:))); ones(size(rho_BA_all(ind_good,:)))];
    XXX2_JC = [rho_JC_all(ind_good,:); rho_JC_all(ind_good,:)]; % rho_JC
    %
    mdl_ABBA = fitlm([XXX1_AB, XXX1_BA, XXX2_JC], YYY_ABBA,'Intercept',false);
    aAB = mdl_ABBA.Coefficients.Estimate(1);
    aBA = mdl_ABBA.Coefficients.Estimate(2);
    bABBA = mdl_ABBA.Coefficients.Estimate(3);
    paAB = coefTest(mdl_ABBA,[1,0,0],0);
    paBA = coefTest(mdl_ABBA,[0,1,0],0);
    pbABBA = coefTest(mdl_ABBA,[0,0,1],1);
    %
    CI95_ABBA = coefCI(mdl_ABBA,0.05);
    CIaAB = CI95_ABBA(1,:);
    CIaABmean = mean(CIaAB);
    CIaABstd = CIaAB(2)-CIaABmean;
    %
    CIaBA = CI95_ABBA(2,:);
    CIaBAmean = mean(CIaBA);
    CIaBAstd = CIaBA(2)-CIaBAmean;
    %
    CIbABBA = CI95_ABBA(3,:);
    CIbABBAmean = mean(CIbABBA);
    CIbABBAstd = CIbABBA(2)-CIbABBAmean;
    %
    subplot(nmonkeys,nplottypes+1,iplottype+1+(imonkey-1)*(nplottypes+1));
    hold on; 
    XXX_JCplot = rho_JC_all(ind_good,:);
    YYY_ABplot = rho_AB_all(ind_good,:);
    YYY_BAplot = rho_BA_all(ind_good,:);
    %
    XXX_lim = [rho_JC_all(ind_good,:); rho_AB_all(ind_good,:); rho_BA_all(ind_good,:)];       
    YYY_lim = XXX_lim;  
%     XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20];
%     YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/20,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/20];
    XX = [0 7.5];
    YY = [0 7.5];
    %    
    YY_AB = XX.*bABBA + aAB;
    hold on; plot(XX,YY_AB,'-','LineWidth',1.5,'color',[1. .2 .2]);
    YY_BA = XX.*bABBA + aBA;
    hold on; plot(XX,YY_BA,'-','LineWidth',1.5,'color',[ .2 .2 1.]);
    %
    hold on; plot(XX,YY,'--','LineWidth',0.5,'color',[0.7 0.7 0.7]);
    %
    hold on; plot(XXX_JCplot, YYY_ABplot, 'o','MarkerSize',6,'color',[1. .2 .2]);
    hold on; plot(XXX_JCplot, YYY_BAplot, 'o','MarkerSize',6,'color',[.2 .2 1.]);
    %
    xlabel('relative value \rho of Task 1');   
    ylabel('relative value \rho of Task 2');  
    legend({['Task 2 AB'], ['Task 2 BA']},'Location','northwest')
    set(gca,'fontsize',11)    
    box off; axis([XX,YY]); axis square  
    %
    h = text(XX(1)+(abs(XX(1))+abs(XX(2)))/2, YY(1)+(abs(YY(1))+abs(YY(2)))/8, ...
                   { ['slope: b2 = '  ,num2str(bABBA,'%.2f'),'±',num2str(CIbABBAstd,'%.2f'), '; p(null=1) = ',  num2str(pbABBA,'%1.1g')];...
                     ['intercept (AB): b1 = ',num2str(aAB,'%.2f'),'±',num2str(CIaABstd,'%.2f'), '; p(null=0) = ',num2str(paAB,'%1.1g')];...
                     ['intercept (BA): b1 = ',num2str(aBA,'%.2f'),'±',num2str(CIaBAstd,'%.2f'), '; p(null=0) = ',num2str(paBA,'%1.1g')];...
                     ['intercept difference: p(ANCOVA) = ',num2str(p_intc,'%1.1g')];
                   },'FontSize',10);
    
    % add monkey names
    if isequal(monkey_ana,monkeys_ana{1})
        axes('position',[.055 .70 .2 .05]);
    elseif isequal(monkey_ana,monkeys_ana{2})
        axes('position',[.055 .25 .2 .05]); 
    end
    hold on; h = text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    set(h,'Rotation',90);
end



% % %
% rho JC AB and BA
% rho ratation angle distribution
if 0 % onlyplotpapers
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1650 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot3512 = length(htt);
    else
        figure(isumplot3512);
    end
    %
    rho_JC_all = rho_all(:,1);
    rho_SO_all = rho_all(:,2);
    rho_AB_all = rhoABBA_all(:,1);
    rho_BA_all = rhoABBA_all(:,2);
    %
    plottypes1 = {'JC','JC'};
    plottypes2 = {'AB','BA'};
    xlabelnames = {'relative value \rho of Task 1',     'relative value \rho of Task 1',   };
    ylabelnames = {'relative value \rho of Task 2(AB)', 'relative value \rho of Task 2(BA)'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes+3,iplottype+(imonkey-1)*(nplottypes+3));
        eval(['XXX = rho_',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = rho_',plottypename2,'_all(ind_good,:);'])
        %
        XXX_lim = [rho_JC_all(ind_good,:); rho_AB_all(ind_good,:); rho_BA_all(ind_good,:)];       
        YYY_lim = XXX_lim;   
        %
        % plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        aa = deming(XXX,YYY);
        % XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5];
        % YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5];
        XX = [0,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20];
        YY = [0,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/20];
        YYfit = aa(2)*XX+aa(1);
        [~, pp_ttest] = ttest(XXX, YYY);
        [pp_wil, ~  ] = signrank(XXX, YYY);
        hold on; plot(XX,YY,'--','LineWidth',1,'color',[0.7 0.7 0.7]);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        plot(XXX, YYY, 'ko','MarkerSize',6,'LineWidth',1);
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/8,...
            {['t test: p = ',num2str(pp_ttest,'%1.1g')]; ...
             ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')]; ...             
             ['mean \Delta = ',num2str(mean(-XXX+YYY),'%.2f')]; ...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions'];...
             }, 'fontsize', 9);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
%         title({[monkey_ana];['analyzed session # =',num2str(sum(ind_good))]});
        % xTickrange = round(XX(1),0):round((round(XX(2),0)-round(XX(1),0))/3):round(XX(2),0);
        xTickrange = round(XX(1),0):2:round(XX(2),0);
        xTickLabelrange = num2str(xTickrange');
        % yTickrange = round(YY(1),0):round((round(YY(2),0)-round(YY(1),0))/3):round(YY(2),0);
        yTickrange = round(YY(1),0):2:round(YY(2),0);
        yTickLabelrange = num2str(yTickrange');
        set(gca, 'XTick', xTickrange, ...                  
                 'XTickLabel', xTickLabelrange, ...                    
                 'YTick', yTickrange, ...
                 'YTickLabel', yTickLabelrange);
        set(gca,'fontsize',11)    
        box off; axis([XX,YY]); axis square  
    end     
    %
    % calculate rotation angle centering at (1,1)
    nvects = length(rho_JC_all);
    thetaJCAB_all = [];
    thetaJCBA_all = [];
    for ivect = 1:nvects
        JCABvect = [rho_JC_all(ivect)-1, rho_AB_all(ivect)-1];
        diagvect = [1,1];
        thetaJCAB = acos(dot(JCABvect,diagvect)/(norm(JCABvect)*norm(diagvect)));
        if rho_JC_all(ivect)>rho_AB_all(ivect), thetaJCAB = -thetaJCAB; end
        if thetaJCAB > pi/2 | thetaJCAB < -pi/2, thetaJCAB = 0; end
        thetaJCAB_all(ivect,1) = thetaJCAB;
        %
        JCBAvect = [rho_JC_all(ivect)-1, rho_BA_all(ivect)-1];
        diagvect = [1,1];
        thetaJCBA = acos(dot(JCBAvect,diagvect)/(norm(JCBAvect)*norm(diagvect)));
        if rho_JC_all(ivect)>rho_BA_all(ivect), thetaJCBA = -thetaJCBA; end
        if thetaJCBA > pi/2 | thetaJCBA < -pi/2, thetaJCBA = 0; end
        thetaJCBA_all(ivect,1) = thetaJCBA;
    end
    %
    plottypes1 = {'thetaJCAB', 'thetaJCBA'};
    xlabelnames = {'\theta (Task 2 (AB) vs Task 1)', '\theta (Task 2 (BA) vs Task 1)'};
    nplottypes2 = length(plottypes1);
    for iplottype = 1:nplottypes2    
        plottypename1 = plottypes1{iplottype};
        xlabelname = xlabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes+nplottypes2+1,iplottype+nplottypes+(imonkey-1)*(nplottypes+nplottypes2+1));
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        %      
        XX = [-pi/2,pi/2];
        if isequal(monkey_ana,'Gervinho')
            YY = [0, 46];
        elseif isequal(monkey_ana,'Juan')
            YY = [0, 18];
        end
        edges = [XX(1):pi/80:XX(2)];
        % histogram(miu_ord_all(ind_good),edges,'FaceColor',[.4 .4 .4]); hold on;
        plot([0 0],YY, '--','Color', [0.5 0.5 0.5],'LineWidth',1); hold on;
        histogram(XXX,edges,'FaceColor',[.4 .4 .4]); hold on;
        [~, pp_ttest] = ttest(XXX);
        [ pp_wil,  ~] = signrank(XXX);
        text(XX(1)+XX(2)/20,YY(2)-YY(2)/8, {['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                             ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                             ['mean = ',num2str(mean(XXX),'%.2f')];...
                                             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']},'fontsize',9);
        xlabel(xlabelname);
        ylabel('N of sessions');
        % title({['monkey ',monkey_ana]});
        axis([XX YY]); axis square; box off
        set(gca,'fontsize',11)   
    end
    %
    % rotation fitting
    YYY = [rho_AB_all(ind_good,:); rho_BA_all(ind_good,:)]-1; % rho_SO-1
    XXX1 = [ones(size(rho_AB_all(ind_good,:))); zeros(size(rho_BA_all(ind_good,:)))];
    XXX2 = [zeros(size(rho_AB_all(ind_good,:)));ones(size(rho_BA_all(ind_good,:)))];
    XXX3 = [rho_JC_all(ind_good,:); rho_JC_all(ind_good,:)]-1;
    XXX = [XXX1, XXX2, XXX3]; % delta_AB, delta_BA, rho_JC-1
    %
    mdl = fitlm(XXX, YYY,'Intercept',false);
    aAB = mdl.Coefficients.Estimate(1);
    aBA = mdl.Coefficients.Estimate(2);
    bbb = mdl.Coefficients.Estimate(3);
    paAB = mdl.Coefficients.pValue(1);
    paBA = mdl.Coefficients.pValue(2);
    pbbb = mdl.Coefficients.pValue(3);
    %
    subplot(nmonkeys,nplottypes+nplottypes2+1,iplottype+nplottypes+1+(imonkey-1)*(nplottypes+nplottypes2+1));
    hold on; 
    h = text(0,0.5,{['\rhoTask2 - 1 = aAB * \deltaAB + aBA * \deltaBA + b * (\rhoTask1 - 1)'];...
                    ['aAB = ',num2str(aAB,'%.2f'),'; p(aAB) = ',num2str(paAB,'%1.1g')];...
                    ['aBA = ',num2str(aBA,'%.2f'),'; p(aBA) = ',num2str(paBA,'%1.1g')];...
                    ['b = '  ,num2str(bbb,'%.2f'),'; p(b) = ',  num2str(pbbb,'%1.1g')]},'FontSize',10);
    axis off
      
    %
    % add monkey names
    if isequal(monkey_ana,'Gervinho')
        axes('position',[.055 .70 .2 .05]);
    elseif isequal(monkey_ana,'Juan')
        axes('position',[.055 .25 .2 .05]); 
    end
    hold on; h = text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    set(h,'Rotation',90);
end



% % %
% rho AB and BA
% steepness AB and BA
if 1
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1650 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot3521 = length(htt);
    else
        figure(isumplot3521);
    end
    rho_AB_all = rho_SOABBA_sep_all(:,1);
    rho_BA_all = rho_SOABBA_sep_all(:,2);
    steep_AB_all = steepness_SOABBA_sep_all(:,1);
    steep_BA_all = steepness_SOABBA_sep_all(:,2);
    plottypes1 = {'rho_AB','steep_AB'};
    plottypes2 = {'rho_BA','steep_BA'};
    xlabelnames = {'relative value \rho of Task 2(AB)', 'steepness \eta of Task 2(AB)'};
    ylabelnames = {'relative value \rho of Task 2(BA)', 'steepness \eta of Task 2(BA)'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes,iplottype+(imonkey-1)*(nplottypes));
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])  
        %
        aa = deming(XXX,YYY);
        XX = [0,(max(XXX))+(abs((min(XXX)))+abs((max(XXX))))/20];
        YY = [0,(max(YYY))+(abs((min(YYY)))+abs((max(YYY))))/20];
        YYfit = aa(2)*XX+aa(1);
        [~, pp_ttest] = ttest(XXX, YYY);
        [pp_wil, ~  ] = signrank(XXX, YYY);
        hold on; plot(XX,YY,'--','LineWidth',1,'color',[0.7 0.7 0.7]);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/8,...
            {['t test: p = ',num2str(pp_ttest,'%1.1g')]; ...
             ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')]; ...             
             ['mean \Delta = ',num2str(mean(-XXX+YYY),'%.2f')]; ...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions'];...
             }, 'fontsize', 11);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        set(gca,'fontsize',13)    
        box off; axis([XX,YY]); axis square  
    end    
    %
    % add monkey names
    if isequal(monkey_ana,'Gervinho')
        axes('position',[.055 .70 .2 .05]);
        hold on; h = text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
        set(h,'Rotation',90);
    elseif isequal(monkey_ana,'Juan')
        axes('position',[.055 .25 .2 .05]); 
        hold on; h = text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
        set(h,'Rotation',90);
    end
end


% % %
% rho JC SO AB and BA _ first half, second half
% if ~onlyplotpapers
if 0
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot32 = length(htt);
    else
        figure(isumplot32);
    end
    rho_JC_all = rho_all(:,1);
    rho_SO_all = rho_all(:,2);
    rho_AB_all = rhoABBA_all(:,1);
    rho_BA_all = rhoABBA_all(:,2);
    ind_half1 = rho_JC_all<=median(rho_JC_all(ind_good,:)) & ind_good;
    ind_half2 = rho_JC_all>median(rho_JC_all(ind_good,:)) & ind_good;
    halfcolor = [0.7,0.1,0.0; 0.0,0.1,0.7];
    plottypes1 = {'JC','JC','JC','AB'};
    plottypes2 = {'SO','AB','BA','BA'};
    xlabelnames = {'\rho Task 1', '\rho Task 1',     '\rho Task 1',     '\rho Task 2(AB)'};
    ylabelnames = {'\rho Task 2', '\rho Task 2(AB)', '\rho Task 2(BA)', '\rho Task 2(BA)'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        eval(['XXX_lim = rho_',plottypename1,'_all(ind_good,:);'])
        eval(['YYY_lim = rho_',plottypename2,'_all(ind_good,:);'])
        XXX_lim = [XXX_lim; YYY_lim];
        YYY_lim = XXX_lim;   
        %
        subplot(nmonkeys,nplottypes,iplottype+(imonkey-1)*nplottypes);
        %
        for ihalf = 1:2    
            eval(['XXX = rho_',plottypename1,'_all(ind_half',num2str(ihalf),',:);'])
            eval(['YYY = rho_',plottypename2,'_all(ind_half',num2str(ihalf),',:);'])   
            %
            plot(XXX, YYY, 'o','MarkerSize',8,'LineWidth',1,'color',halfcolor(ihalf,:));
            eval(['aa',num2str(ihalf),' = deming(XXX,YYY);'])
            XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5];
            YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5];
            YYfit = aa(2)*XX+aa(1);
            eval(['[~, pp_ttest',num2str(ihalf),'] = ttest(XXX, YYY);'])
            eval(['[pp_wil',num2str(ihalf),', ~  ] = signrank(XXX, YYY);'])
            hold on; plot(XX,YY,'--','LineWidth',1,'color',[0.7 0.7 0.7]);
            Sigma_ell = cov(XXX, YYY);
            mu_ell(1) = mean(XXX);
            mu_ell(2) = mean(YYY);  
            hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);   
        end
        text(XX(1)+(abs((min(XXX)))+abs((max(XXX))))/12,YY(2)-(abs((min(YYY)))+abs((max(YYY))))/8,...
            {['low half - t test: p = ',num2str(pp_ttest1,'%1.1g')]; ...
             ['low half - Wilcoxon: p = ',num2str(pp_wil1,'%1.1g')]; ...
             ['low half - deming slopes = ',num2str(aa1(2),'%.4f')]; ...
             ['high half - t test: p = ',num2str(pp_ttest2,'%1.1g')]; ...
             ['high half - Wilcoxon: p = ',num2str(pp_wil2,'%1.1g')]; ...
             ['high half - deming slopes = ',num2str(aa2(2),'%.4f')]; ...
             }, 'fontsize', 8);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        title({[monkey_ana];['analyzed session # =',num2str(sum(ind_good))]});
        box off; axis([XX,YY]); axis square  
    end
end


% % % 
% steepness difference and relative value difference
% steepness difference and miu
% relative value different and miu
if onlyplotpapers
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot5 = length(htt);
    else
        figure(isumplot5);
    end
    % rho_delta_all = rho_all(:,2) - rho_all(:,1);
    rho_delta_all = (rho_all(:,2) - rho_all(:,1))./(rho_all(:,2) + rho_all(:,1));
    steep_delta_all = steepness_all(:,2) - steepness_all(:,1);
    miu_rho_all = 2.*miu_ord_all.*rho_SO_all;
    plottypes1 = {'miu_rho','miu_rho','rho_delta'};
    plottypes2 = {'rho_delta','steep_delta','steep_delta'};
    xlabelnames = {'order bias \epsilon','order bias \epsilon', 'normalized \Delta\rho (Task 2 - Task 1)'};
    ylabelnames = {'normalized \Delta\rho (Task 2 - Task 1)','\Delta\eta (Task 2 - Task 1)', '\Delta\eta (Task 2 - Task 1)'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes,iplottype+(imonkey-1)*nplottypes);
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])
        %
        eval(['XXX_lim = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY_lim = ',plottypename2,'_all(ind_good,:);'])
        %
        % plot(XXX, YYY, 'ko','MarkerSize',10,'LineWidth',1);
        aa = deming(XXX,YYY);
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',4,'color',[0.7 0.7 0.7]);
        hold on; plot(XX,[0 0],'--','LineWidth',1,'color',[0.2 0.2 0.2]);
        hold on; plot([0 0],YY,'--','LineWidth',1,'color',[0.2 0.2 0.2]);
        plot(XXX, YYY, 'ko','MarkerSize',10,'LineWidth',1);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/10,YY(2)-(abs((min(YY)))+abs((max(YY))))/20,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];  ...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions'];...
             }, 'fontsize', 11);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);        
        box off; axis([XX,YY]); axis square  
        set(gca,'fontsize',13)
    end    
    %
    % add monkey names
    if isequal(monkey_ana,monkeys_ana{1})
        axes('position',[.055 .70 .2 .05]);
    elseif isequal(monkey_ana,monkeys_ana{2})
        axes('position',[.055 .25 .2 .05]); 
    end
    hold on; h = text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    set(h,'Rotation',90);
end


% order bias v.s. steepness difference
% rho SO v.s. steepness difference
% rho SO v.s. order bias
if 0
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1650 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot53 = length(htt);
    else
        figure(isumplot53);
    end
    %
    rho_JC_all = rho_all(:,1);
    rho_SO_all = rho_all(:,2);
    miu_rho_all = 2.*miu_ord_all.*rho_SO_all;
    rho_delta_all = rho_all(:,2) - rho_all(:,1);
    steep_delta_all = steepness_all(:,2) - steepness_all(:,1);
    plottypes1 = {'miu_rho',    'rho_SO',      'rho_SO'};
    plottypes2 = {'steep_delta','steep_delta', 'miu_rho'};
    xlabelnames = { 'order bias \epsilon of Task 2',          'relative value \rho of Task 2',           'relative value \rho of Task 2'};
    ylabelnames = { '\Delta steepness \eta (Task 2 - Task 1)','\Delta steepness \eta (Task 2 - Task 1)', 'order bias \epsilon of Task 2'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes,iplottype+(imonkey-1)*nplottypes);
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])
        %
        eval(['XXX_lim = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY_lim = ',plottypename2,'_all(ind_good,:);'])
        %
        % plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        aa = deming(XXX,YYY);
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/2];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',4,'color',[0.7 0.7 0.7]);
        plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/8,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']}, 'fontsize', 11);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        set(gca,'fontsize',13)    
        box off; axis([XX,YY]); axis square  
    end    
    %
    % add monkey names
    if isequal(monkey_ana,'Gervinho')
        axes('position',[.055 .70 .2 .05]);
    elseif isequal(monkey_ana,'Juan')
        axes('position',[.055 .25 .2 .05]); 
    end
    hold on; h = text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    set(h,'Rotation',90);
end


% % %
% behavioral effects v.s. fraction of choosing A
if dofractofA
    JCSOs = {'JC','SO'};
    Task12s = {'Task 1', 'Task 2'};
    for iJCSO = 1:2
    JCSO = JCSOs{iJCSO};
    Task12 = Task12s{iJCSO};
    %
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        eval(['isumplot',num2str(678+iJCSO),' = length(htt);'])
    else
        eval(['figure(isumplot',num2str(678+iJCSO),');'])
    end
    %
    rho_JC_all = rho_all(:,1);
    rho_SO_all = rho_all(:,2);
    fractofA_JC_all = fractofA_all(:,1);
    fractofA_SO_all = fractofA_all(:,2);
    steep_JC_all = steepness_all(:,1);
    steep_SO_all = steepness_all(:,2);
    % rho_delta_all = rho_all(:,2) - rho_all(:,1);
    rho_delta_all = (rho_all(:,2) - rho_all(:,1))./(rho_all(:,2) + rho_all(:,1));
    steep_delta_all = steepness_all(:,2) - steepness_all(:,1);
    miu_rho_all = 2.*miu_ord_all.*rho_SO_all;
    eval(['plottypes1 = repelem({''fractofA_',JCSO,'''},7,1);'])
    plottypes2 = {'steep_JC', 'steep_SO', 'rho_JC', 'rho_SO', 'miu_rho', 'steep_delta', 'rho_delta'};
    eval(['xlabelnames = repelem({''fraction of choosing A in ',Task12,'''},7,1);'])
    ylabelnames = {'steepness \eta of Task 1', 'steepness \eta of Task 2', 'relative value \rho of Task 1', 'relative value \rho of Task 2', ...
                   'order bias \epsilon of Task 2', '\Delta\eta (Task 2 - Task 1)', 'normalized \Delta\rho (Task 2 - Task 1)'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys*ceil(nplottypes/5),5,iplottype+(imonkey-1)*ceil(nplottypes/5)*5);
        % subplot(nmonkeys,nplottypes,iplottype+(imonkey-1)*nplottypes);
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])
        %
        eval(['XXX_lim = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY_lim = ',plottypename2,'_all(ind_good,:);'])
        %
        % plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        aa = deming(XXX,YYY);
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/2];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',2,'color',[0.7 0.7 0.7]);
        plot(XXX, YYY, 'ko','MarkerSize',5,'LineWidth',1);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/5,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']}, 'fontsize', 8);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        set(gca,'fontsize',9)    
        box off; axis([XX,YY]); axis square  
    end    
    %
    % add monkey names
    if isequal(monkey_ana,'Gervinho')
        axes('position',[.025 .25 .2 .05]);
    elseif isequal(monkey_ana,'Juan')
        axes('position',[.025 .75 .2 .05]); 
    end
    text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    end
end

% % %
% behavioral effects v.s. fraction of choosing A
if domeanvalues
    JCSOs = {'JC','SO'};
    Task12s = {'Task 1', 'Task 2'};
    for iJCSO = 1:2
    JCSO = JCSOs{iJCSO};
    Task12 = Task12s{iJCSO};
    %
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        eval(['isumplot',num2str(788+iJCSO),' = length(htt);'])
    else
        eval(['figure(isumplot',num2str(788+iJCSO),');'])
    end
    %
    rho_JC_all = rho_all(:,1);
    rho_SO_all = rho_all(:,2);
    fractofA_JC_all = fractofA_all(:,1);
    fractofA_SO_all = fractofA_all(:,2);
    meanOVA_JC_all = meanOVA_all(:,1);
    meanOVB_JC_all = meanOVB_all(:,1);
    meanOVA_SO_all = meanOVA_all(:,2);
    meanOVB_SO_all = meanOVB_all(:,2);
    meanCVA_JC_all = meanCVA_all(:,1);
    meanCVB_JC_all = meanCVB_all(:,1);
    meanCVA_SO_all = meanCVA_all(:,2);
    meanCVB_SO_all = meanCVB_all(:,2);
    deltaOV_JC_all = meanOVA_JC_all - meanOVB_JC_all;
    deltaOV_SO_all = meanOVA_SO_all - meanOVB_SO_all;
    deltaCV_JC_all = meanCVA_JC_all - meanCVB_JC_all;
    deltaCV_SO_all = meanCVA_SO_all - meanCVB_SO_all;
    rho_delta_all = 2*(rho_all(:,2) - rho_all(:,1))./(rho_all(:,2) + rho_all(:,1));
    plottypes1 = {['fractofA_',JCSO], ['deltaOV_',JCSO], ['deltaCV_',JCSO]};
    plottypes2 = {'rho_delta', 'rho_delta', 'rho_delta'};
    eval(['xlabelnames = repelem({''fraction of choosing A in ',Task12,'''},1,1);'])
    xlabelnames = {['fraction of choosing A in ',Task12], ['<OVA> - <OVB> in ',Task12], ['<CVA> - <CVB> in ',Task12]};
    ylabelnames = {'Preference Bias Index', 'Preference Bias Index', 'Preference Bias Index'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot((nmonkeys+2)*ceil(nplottypes/3),3,iplottype+(imonkey-1)*2*ceil(nplottypes/3)*3);
        % subplot(nmonkeys,nplottypes,iplottype+(imonkey-1)*nplottypes);
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])       
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])
%         ind_outlier = XXX>12;
%         XXX(ind_outlier) = [];
%         YYY(ind_outlier) = [];
        %
        XXX_lim = XXX;
        YYY_lim = YYY;
        %
        % plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        aa = deming(XXX,YYY);
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/2];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',1.5,'color',[0.7 0.7 0.7]);
        hold on; plot(XX,[0 0],'--','LineWidth',0.5,'color',[0.0 0.0 0.0]);
        if iplottype == 1
            hold on; plot([.5 .5],YY,'--','LineWidth',0.5,'color',[0.0 0.0 0.0]);
        else
            hold on; plot([0 0],YY,'--','LineWidth',0.5,'color',[0.0 0.0 0.0]);
        end
        plot(XXX, YYY, 'ko','MarkerSize',5,'LineWidth',1);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/5,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']}, 'fontsize', 8);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        set(gca,'fontsize',9)    
        box off; axis([XX,YY]); axis square  
        % % % 
        % % %
        %
        subplot((nmonkeys+2)*ceil(nplottypes/3),3,iplottype+((imonkey-1)*2+1)*ceil(nplottypes/3)*3);
        edges = [XX(1):(XX(2)-XX(1))/20:XX(2)];
        HHH = histogram(XXX,edges,'FaceColor',[.4 .4 .4]); hold on;  
        YY = [0 max(HHH.Values)*1.25];
        if iplottype == 1
            [~, pp_ttest] = ttest(XXX-0.5);
            [ pp_wil,  ~] = signrank(XXX-0.5);
            plot([0.5 0.5],YY, '--','Color', [0.0 0.0 0.0],'LineWidth',0.5); hold on;
        else
            [~, pp_ttest] = ttest(XXX);
            [ pp_wil,  ~] = signrank(XXX);
            plot([0 0],YY, '--','Color', [0.0 0.0 0.0],'LineWidth',0.5); hold on;
        end
        text(XX(1)+XX(2)/1.5,YY(2)-YY(2)/8, {['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                             ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                             ['mean = ',num2str(mean(XXX),'%.2f')];...
                                             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']},'fontsize',8);
        xlabel([xlabelname]);   
        ylabel('N of sessions');
        axis([XX YY]); axis square; box off
        set(gca,'fontsize',9)   

    end    
    %
    % add monkey names
    if isequal(monkey_ana,'Gervinho')
        axes('position',[.025 .25 .2 .05]);
    elseif isequal(monkey_ana,'Juan')
        axes('position',[.025 .75 .2 .05]); 
    end
    text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    end
end

% % %
% behavioral effects v.s. session value ranges 
if dovaluerange
    JCSOs = {'JC','SO'};
    Task12s = {'Task 1', 'Task 2'};
    for iJCSO = 1:2
    JCSO = JCSOs{iJCSO};
    Task12 = Task12s{iJCSO};
    %
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        eval(['isumplot',num2str(555+iJCSO),' = length(htt);'])
    else
        eval(['figure(isumplot',num2str(555+iJCSO),');'])
    end
    %
    rho_JC_all = rho_all(:,1);
    rho_SO_all = rho_all(:,2);
    eval(['sessionrange_all = sqrt(rho_',JCSO,'_all.*quantrange_all(:,1).*quantrange_all(:,2));'])    
    steep_JC_all = steepness_all(:,1);
    steep_SO_all = steepness_all(:,2);
    rho_delta_all = rho_all(:,2) - rho_all(:,1);
    steep_delta_all = steepness_all(:,2) - steepness_all(:,1);
    miu_rho_all = 2.*miu_ord_all.*rho_SO_all;
    plottypes1 = repelem({'sessionrange'},7,1);
    plottypes2 = {'steep_JC', 'steep_SO', 'rho_JC', 'rho_SO', 'miu_rho', 'steep_delta', 'rho_delta'};
    eval(['xlabelnames = repelem({''session range in ',Task12,'''},7,1);'])
    ylabelnames = {'steepness \eta of Task 1', 'steepness \eta of Task 2', 'relative value \rho of Task 1', 'relative value \rho of Task 2',...
                   'order bias \epsilon of Task 2','\Delta steepness \eta (Task 2 - Task 1)', '\Delta relative value \rho (Task 2 - Task 1)'};
    nplottypes = length(plottypes1);   
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys*ceil(nplottypes/5),5,iplottype+(imonkey-1)*ceil(nplottypes/5)*5);
        % subplot(nmonkeys,nplottypes,iplottype+(imonkey-1)*nplottypes);
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])
        %
        eval(['XXX_lim = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY_lim = ',plottypename2,'_all(ind_good,:);'])
        %
        % plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        aa = deming(XXX,YYY);
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/2];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',2,'color',[0.7 0.7 0.7]);
        plot(XXX, YYY, 'ko','MarkerSize',5,'LineWidth',1);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/5,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']}, 'fontsize', 8);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        set(gca,'fontsize',9)    
        box off; axis([XX,YY]); axis square  
    end    
    %
    % add monkey names
    if isequal(monkey_ana,'Gervinho')
        axes('position',[.025 .75 .2 .05]);
    elseif isequal(monkey_ana,'Juan')
        axes('position',[.025 .25 .2 .05]); 
    end
    text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    end
end


% % %
% behavioral effects v.s. value ranges difference 
if 0 % dovaluerange
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot551 = length(htt);
    else
        figure(isumplot551);
    end
    %
    rho_JC_all = rho_all(:,1);
    rangediff_all = rho_JC_all.*quantrange_all(:,1) - quantrange_all(:,2);    
    steep_JC_all = steepness_all(:,1);
    rho_delta_all = rho_all(:,2) - rho_all(:,1);
    steep_delta_all = steepness_all(:,2) - steepness_all(:,1);
    plottypes1 = {'rangediff', 'rangediff', 'rangediff', 'rangediff',   'rangediff'};
    plottypes2 = {'steep_JC',  'rho_JC',    'miu_ord',   'steep_delta', 'rho_delta'};
    xlabelnames = {'\DeltaVA - \DeltaVB',      '\DeltaVA - \DeltaVB',           '\DeltaVA - \DeltaVB',           '\DeltaVA - \DeltaVB',                     '\DeltaVA - \DeltaVB'};
    ylabelnames = {'steepness \eta of Task 1', 'relative value \rho of Task 1', 'order bias \epsilon of Task 2', '\Delta steepness \eta (Task 2 - Task 1)', '\Delta relative value \rho (Task 2 - Task 1)'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes,iplottype+(imonkey-1)*nplottypes);
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])
        %
        eval(['XXX_lim = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY_lim = ',plottypename2,'_all(ind_good,:);'])
        %
        % plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        aa = deming(XXX,YYY);
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/2];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',4,'color',[0.7 0.7 0.7]);
        plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/5,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']}, 'fontsize', 11);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        set(gca,'fontsize',13)    
        box off; axis([XX,YY]); axis square  
    end    
    %
    % add monkey names
    if isequal(monkey_ana,monkeys_ana{1})
        axes('position',[.025 .75 .2 .05]);
    elseif isequal(monkey_ana,monkeys_ana{2})
        axes('position',[.025 .25 .2 .05]); 
    end
    text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
end


% % %
% plot rho effects separate for different juice pair
if plotjuicepairs
    %
    nsessions = length(juicetypes_all);
    juicetypes_merge = {};
    for isession = 1:nsessions
        juicetypes_merge(isession,:) = {[juicetypes_all{isession,1}, juicetypes_all{isession,2}]};
    end
    juicepair_uni_merge = unique(juicetypes_merge);
    juicepair_uni = {};
    juicepair_lgt = {};
    njuicepair_uni = length(juicepair_uni_merge);
    for ipair_uni = 1:njuicepair_uni
        juicepair_merge = juicepair_uni_merge{ipair_uni};
        ind_pair = ismember(juicetypes_merge,juicepair_merge);
        juicepair_unmarge = juicetypes_all(ind_pair,:);
        juicepair_unmarge = juicepair_unmarge(1,:);
        juicepair_uni(ipair_uni,:) = juicepair_unmarge;
        juicepair_lgt(ipair_uni,:) = {[juicepair_unmarge{1}, '; ' ,juicepair_unmarge{2}]};
        rho_uni(ipair_uni,:) = nanmean(rho_JC_all(ind_pair,:));
    end
    [~,juicepair_index] = ismember(juicetypes_merge,juicepair_uni_merge);  
    juicepair_index = juicepair_index(ind_good,:);
    %
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1650 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot352 = length(htt);
    else
        figure(isumplot352);
    end
    rho_JC_all = rho_all(:,1);
    rho_SO_all = rho_all(:,2);
    rho_AB_all = rhoABBA_all(:,1);
    rho_BA_all = rhoABBA_all(:,2);
    plottypes1 = {'JC','JC'};
    plottypes2 = {'AB','BA'};
    xlabelnames = {'relative value \rho of Task 1',     'relative value \rho of Task 1',   };
    ylabelnames = {'relative value \rho of Task 2(AB)', 'relative value \rho of Task 2(BA)'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes+1,iplottype+(imonkey-1)*(nplottypes+1));
        eval(['XXX = rho_',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = rho_',plottypename2,'_all(ind_good,:);'])
        %
        for ijuicepair = 1:njuicepair_uni
            ind_ipair = juicepair_index == ijuicepair;
            hold on
            plot(XXX(ind_ipair), YYY(ind_ipair), 'o','MarkerSize',7,'LineWidth',1);
        end
        %
        XXX_lim = [rho_JC_all(ind_good,:); rho_AB_all(ind_good,:); rho_BA_all(ind_good,:)];       
        YYY_lim = XXX_lim;   
        aa = deming(XXX,YYY);
        XX = [0,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20];
        YY = [0,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/20];
        YYfit = aa(2)*XX+aa(1);
        [~, pp_ttest] = ttest(XXX, YYY);
        [pp_wil, ~  ] = signrank(XXX, YYY);
        hold on; plot(XX,YY,'--','LineWidth',1,'color',[0.7 0.7 0.7]);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);  
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/8,...
            {['t test: p = ',num2str(pp_ttest,'%1.1g')]; ...
             ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')]; ...             
             ['mean \Delta = ',num2str(mean(-XXX+YYY),'%.2f')]; ...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions'];...
             }, 'fontsize', 11);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        xTickrange = round(XX(1),0):2:round(XX(2),0);
        xTickLabelrange = num2str(xTickrange');
        yTickrange = round(YY(1),0):2:round(YY(2),0);
        yTickLabelrange = num2str(yTickrange');
        set(gca, 'XTick', xTickrange, ...                  
                 'XTickLabel', xTickLabelrange, ...                    
                 'YTick', yTickrange, ...
                 'YTickLabel', yTickLabelrange);
        set(gca,'fontsize',13)    
        box off; axis([XX,YY]); axis square  
    end    
    %
    % plot rho v.s. rho difference
    rho_JC_all = rho_all(:,1);
    rho_delta_all = rho_all(:,2) - rho_all(:,1);
    steep_delta_all = steepness_all(:,2) - steepness_all(:,1);
    plottypes1 = {'rho_JC',};
    plottypes2 = {'rho_delta'};
    xlabelnames = {'relative value \rho of Task 1'};
    ylabelnames = {'\Delta relative value \rho (Task 2 - Task 1)'};
    nplottypes2 = length(plottypes1);
    for iplottype = 1:nplottypes2    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes+1,iplottype+nplottypes+(imonkey-1)*(nplottypes+1));
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])
        %
        for ijuicepair = 1:njuicepair_uni
            ind_ipair = juicepair_index == ijuicepair;
            hold on 
            plot(XXX(ind_ipair), YYY(ind_ipair), 'o','MarkerSize',7,'LineWidth',1);
        end
        %
        eval(['XXX_lim = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY_lim = ',plottypename2,'_all(ind_good,:);'])
        %
        aa = deming(XXX,YYY);
        XX = [0,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/20];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/2];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',4,'color',[0.7 0.7 0.7]);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/8,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']}, 'fontsize', 11);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        xTickrange = round(XX(1),1):2:round(XX(2),1);
        xTickLabelrange = num2str(xTickrange');
        set(gca, 'XTick', xTickrange, ...                  
                 'XTickLabel', xTickLabelrange);
        set(gca,'fontsize',13)    
        box off; axis([XX,YY]); axis square  
        lgt = legend(juicepair_lgt);
        lgt.FontSize = 6;
    end      
    %
    % add monkey names
    if isequal(monkey_ana,'Gervinho')
        axes('position',[.025 .75 .2 .05]);
    elseif isequal(monkey_ana,'Juan')
        axes('position',[.025 .25 .2 .05]); 
    end
    text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
end


% % % 
%plot for reaction time: JC AB BA SO
if doreacttime & ~onlyplotpapers
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot6 = length(htt);
    else
        figure(isumplot6);
    end
    % reaction time
    RT_anatype1 = {'JC','AB','JC','JC'};
    RT_anatype2 = {'SO','BA','AB','BA'};
    xlabelnames = {'RT Task 1', 'RT Task 2(AB)', 'RT Task 1',     'RT Task 1'};
    ylabelnames = {'RT Task 2', 'RT Task 2(BA)', 'RT Task 2(AB)', 'RT Task 2(BA)'};
    nplots = length(RT_anatype1);
    for iplot = 1:nplots   
        subplot(nmonkeys,nplots,iplot+(imonkey-1)*nplots);
        eval(['XXX = RT_',RT_anatype1{iplot},'_all(ind_good,1);'])
        eval(['YYY = RT_',RT_anatype2{iplot},'_all(ind_good,1);'])
        xlabelname = xlabelnames{iplot};
        ylabelname = ylabelnames{iplot};
        plot(XXX,YYY, 'ko','MarkerSize',10);
        hold on 
        minX = floor(min([min(XXX),min(YYY)]))-10;
        maxX = ceil(max([max(XXX),max(YYY)]))+10;
        plot([minX,maxX],[minX,maxX],'--','LineWidth',1,'color',[0.7 0.7 0.7]);
        Sigma_ell = cov(XXX,YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
        [~, pp_ttest] = ttest(XXX,YYY);
        [pp_wil, ~  ] = signrank(XXX,YYY);
        text(minX+(abs(minX)+abs(maxX))/100,maxX-(abs(minX)+abs(maxX))/100,...
             {['t test: p = ',num2str(pp_ttest,'%1.1g')];['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')]},'FontSize',8);
        xlabel(xlabelname);
        ylabel(ylabelname);
        title({['reaction time; Monkey ',monkey_ana];['analyzed session # =',num2str(sum(ind_good))]});
        box off; axis([minX,maxX,minX,maxX]); axis square
    end
end


% % % 
% reaction time and relative value and session value
if doreacttime & dovaluerange
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1650 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot611 = length(htt);
    else
        figure(isumplot611);
    end
    %
    %
    rho_JC_all = rho_all(:,1);
    sessionrange_all = sqrt(rho_JC_all.*quantrange_all(:,1).*quantrange_all(:,2));   
    steep_JC_all = steepness_all(:,1);
    rho_delta_all = rho_all(:,2) - rho_all(:,1);
    steep_delta_all = steepness_all(:,2) - steepness_all(:,1);
    plottypes1 = {'rho_JC', 'sessionrange', 'rho_SO', 'sessionrange'};
    plottypes2 = {'RT_JC',  'RT_JC',        'RT_SO',  'RT_SO'};
    xlabelnames = {'relative value \rho of Task 1', 'session value range', 'relative value \rho of Task 2', 'session value range'};
    ylabelnames = {'RT of Task 1',                  'RT of Task 1',        'RT of Task 2',                  'RT of Task 2'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(nmonkeys,nplottypes,iplottype+(imonkey-1)*nplottypes);
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])
        %
        eval(['XXX_lim = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY_lim = ',plottypename2,'_all(ind_good,:);'])
        %
        % plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        aa = deming(XXX,YYY);
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/3];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/30,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/30];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',4,'color',[0.7 0.7 0.7]);
        plot(XXX, YYY, 'ko','MarkerSize',9,'LineWidth',1);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XX)))+abs((max(XX))))/50,YY(2)-(abs((min(YY)))+abs((max(YY))))/50,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];...
             ['N = ',num2str(nsessions_ana-nsession_steepnessout),' sessions']}, 'fontsize', 11);
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        set(gca,'fontsize',13)    
        box off; axis([XX,YY]); axis square  
    end    
    %
    % add monkey names
    if isequal(monkey_ana,'Gervinho')
        axes('position',[.025 .75 .2 .05]);
    elseif isequal(monkey_ana,'Juan')
        axes('position',[.025 .25 .2 .05]); 
    end
    text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
end


% % % 
% detailed raction time analysis: chA vs chB in JC AB BA SO
if 0 % doreacttime & ~onlyplotpapers
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot62 = length(htt);
    else
        figure(isumplot62);
    end
    % reaction time
    RT_anatype1 = {'JCchA','SOchA','ABchA','BAchA'};
    RT_anatype2 = {'JCchB','SOchB','ABchB','BAchB'};
    xlabelnames = {'RT Task 1(A chosen)', 'RT Task 2(A chosen)', 'RT Task 2(AB; A chosen)', 'RT Task 2(BA; A chosen)'};
    ylabelnames = {'RT Task 1(B chosen)', 'RT Task 2(B chosen)', 'RT Task 2(AB; B chosen)', 'RT Task 2(BA; B chosen)'};
    nplots = length(RT_anatype1);
    for iplot = 1:nplots   
        subplot(nmonkeys,nplots,iplot+(imonkey-1)*nplots);
        eval(['XXX = RT_',RT_anatype1{iplot},'_all(ind_good,1);'])
        eval(['YYY = RT_',RT_anatype2{iplot},'_all(ind_good,1);'])
        xlabelname = xlabelnames{iplot};
        ylabelname = ylabelnames{iplot};
        plot(XXX,YYY, 'ko','MarkerSize',10);
        hold on 
        minX = floor(min([min(XXX),min(YYY)]))-10;
        maxX = ceil(max([max(XXX),max(YYY)]))+10;
        plot([minX,maxX],[minX,maxX],'--','LineWidth',1,'color',[0.7 0.7 0.7]);
        Sigma_ell = cov(XXX,YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);
        hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
        [~, pp_ttest] = ttest(XXX,YYY);
        [pp_wil, ~  ] = signrank(XXX,YYY);
        text(minX+(abs(minX)+abs(maxX))/100,maxX-(abs(minX)+abs(maxX))/100,...
             {['t test: p = ',num2str(pp_ttest,'%1.1g')];['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
              ['mean \Delta = ',num2str(mean(-XXX+YYY),'%.2f')]},'FontSize',8);
        xlabel(xlabelname);
        ylabel(ylabelname);
        title({['reaction time; Monkey ',monkey_ana];['analyzed session # =',num2str(sum(ind_good))]});
        box off; axis([minX,maxX,minX,maxX]); axis square
    end
end


% % % 
%reaction time difference and relative value different and steepness difference
if doreacttime & ~onlyplotpapers
    if imonkey == 1
        figure;
        set(gcf,'position',[110 65 1850 1250], 'PaperPositionMode','auto')
        htt =  findobj('type','figure');
        isumplot7 = length(htt);
    else
        figure(isumplot7);
    end
    RT_all = [RT_JC_all, RT_SO_all];
    % reaction time
    anatype1 = {'RT',  'RT',        'rho'      };
    anatype2 = {'rho', 'steepness', 'steepness'};
    nplots = length(anatype1);
    for iplot = 1:nplots   
        subplot(nmonkeys,nplots,iplot+(imonkey-1)*nplots);
        eval(['XXX = ',anatype1{iplot},'_all(ind_good,1)-',anatype1{iplot},'_all(ind_good,2);'])
        eval(['YYY = ',anatype2{iplot},'_all(ind_good,1)-',anatype2{iplot},'_all(ind_good,2);'])
        % plot(XXX,YYY, 'ko','MarkerSize',10);
        hold on 
        aa = deming(XXX,YYY);
        XX = [floor(min(XXX))-(abs(floor(min(XXX)))+abs(ceil(max(XXX))))/5,ceil(max(XXX))+(abs(floor(min(XXX)))+abs(ceil(max(XXX))))/5];
        YY = [floor(min(YYY))-(abs(floor(min(YYY)))+abs(ceil(max(YYY))))/5,ceil(max(YYY))+(abs(floor(min(YYY)))+abs(ceil(max(YYY))))/5];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        plot(XX,YYfit,'-','LineWidth',4,'color',[0.7 0.7 0.7]);
        plot(XXX,YYY, 'ko','MarkerSize',10);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);    
        text(XX(1)+(abs((min(XXX)))+abs((max(XXX))))/10,YY(2)-(abs((min(YYY)))+abs((max(YYY))))/10,...
            {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
             ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')]}, 'fontsize', 9);
        xlabel(['\Delta',anatype1{iplot},'(Task 1 - Task 2)']);
        ylabel(['\Delta',anatype2{iplot},'(Task 1 - Task 2)']);
        title({[monkey_ana];['analyzed session # =',num2str(sum(ind_good))]});
        box off; axis([XX,YY]); axis square
    end
end


% % %
% choice hysteresis (SO for SO and JC for JC )
if 1
    figure;
    set(gcf,'position',[110 65 1650 1250], 'PaperPositionMode','auto')
    hyst_new_all = 2.*hyst_all.*rho_all;
%     hyst_new_all = hyst_all;
    %
    subplot(2,3,1);
    XX = [-2, 3.5]; 
    if isequal(monkey_ana,'Gervinho')
        YY = [0, 38];
    elseif isequal(monkey_ana,'Juan')
        YY = [0, 25];
    end
    edges = [XX(1):0.1:XX(2)];
    % histogram(miu_ord_all(ind_good),edges,'FaceColor',[.4 .4 .4]); hold on;
    plot([0 0],YY, '--','Color', [0.5 0.5 0.5],'LineWidth',1); hold on;
    histogram(hyst_new_all(ind_good,1),edges,'FaceColor',[.4 .4 .4]); hold on;
    [~, pp_ttest] = ttest(hyst_new_all(ind_good,1));
    [ pp_wil,  ~] = signrank(hyst_new_all(ind_good,1));
    text(XX(1)+XX(2)/1.5,YY(2)-YY(2)/8, {['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                         ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                         ['mean = ',num2str(mean(hyst_new_all(ind_good,1)),'%.2f')];...
                                         ['N = ',num2str(sum(ind_good)),' sessions']},'fontsize',10);
    xlabel({['choice hysteresis']});
    ylabel('N of sessions');
    title(['Task 1 on Task 1']);
    axis([XX YY]); axis square; box off
    set(gca,'fontsize',13)   
    %
    subplot(2,3,2);
    plot([0 0],YY, '--','Color', [0.5 0.5 0.5],'LineWidth',1); hold on;
    histogram(hyst_new_all(ind_good,2),edges,'FaceColor',[.4 .4 .4]); hold on;
    [~, pp_ttest] = ttest(hyst_new_all(ind_good,2));
    [ pp_wil,  ~] = signrank(hyst_new_all(ind_good,2));
    text(XX(1)+XX(2)/1.5,YY(2)-YY(2)/8, {['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                         ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                         ['mean = ',num2str(mean(hyst_new_all(ind_good,2)),'%.2f')];...
                                         ['N = ',num2str(sum(ind_good)),' sessions']},'fontsize',10);
    xlabel({['choice hysteresis']});
    ylabel('N of sessions');
    title(['Task 2 on Task 2']);
    axis([XX YY]); axis square; box off
    set(gca,'fontsize',13) 
    
    % hysteresis and other bhv effects
    %
    % rho_delta_all = rho_all(:,2) - rho_all(:,1);
    rho_delta_all = (rho_all(:,2) - rho_all(:,1))./(rho_all(:,2) + rho_all(:,1));
    steep_delta_all = steepness_all(:,2) - steepness_all(:,1);
    steep_SO_all = steepness_all(:,1);
    miu_rho_all = 2.*miu_ord_all.*rho_all(:,2);
    hyst_delta_all = hyst_new_all(:,2) - hyst_new_all(:,1);
    hyst_JC_all = hyst_new_all(:,1);
    hyst_SO_all = hyst_new_all(:,2);
    plottypes1 = {'hyst_JC', 'miu_rho', 'steep_delta', 'rho_delta'};
%     plottypes2 = {'hyst_SO', 'hyst_SO', 'hyst_SO',     'hyst_SO'};
    plottypes2 = {'hyst_SO', 'hyst_delta', 'hyst_delta',     'hyst_delta'};
    % plottypes2 = {'hyst_SO', 'hyst_JC', 'hyst_JC',     'hyst_JC'};
    xlabelnames = {'choice hysteresis of Task 1', 'order bias \epsilon',         '\Delta\eta (Task 2 - Task 1)', 'normalized \Delta\rho (Task 2 - Task 1)'};
%     ylabelnames = {'choice hysteresis of Task 2', 'choice hysteresis of Task 2', 'choice hysteresis of Task 2',  'choice hysteresis of Task 2'};
    ylabelnames = {'choice hysteresis of Task 2', '\Delta\xi (Task 2 - Task 1)', '\Delta\xi (Task 2 - Task 1)',  '\Delta\xi (Task 2 - Task 1)'};
    % ylabelnames = {'choice hysteresis of Task 2', 'choice hysteresis of Task 1', 'choice hysteresis of Task 1',  'choice hysteresis of Task 1'};
    nplottypes = length(plottypes1);
    for iplottype = 1:nplottypes    
        plottypename1 = plottypes1{iplottype};
        plottypename2 = plottypes2{iplottype};
        xlabelname = xlabelnames{iplottype};
        ylabelname = ylabelnames{iplottype};
        %
        subplot(2,3,iplottype+2);
        eval(['XXX = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY = ',plottypename2,'_all(ind_good,:);'])
        %
        eval(['XXX_lim = ',plottypename1,'_all(ind_good,:);'])
        eval(['YYY_lim = ',plottypename2,'_all(ind_good,:);'])
        %
        % plot(XXX, YYY, 'ko','MarkerSize',10,'LineWidth',1);
        aa = deming(XXX,YYY);
        XX = [(min(XXX_lim))-(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5,(max(XXX_lim))+(abs((min(XXX_lim)))+abs((max(XXX_lim))))/5];
        YY = [(min(YYY_lim))-(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5,(max(YYY_lim))+(abs((min(YYY_lim)))+abs((max(YYY_lim))))/5];
        YYfit = aa(2)*XX+aa(1);
        [RR_Spe,pp_Spe] = corr(XXX,YYY,'Type','Spearman');
        [RR_Pea,pp_Pea] = corr(XXX,YYY,'Type','Pearson');
        hold on; plot(XX,YYfit,'-','LineWidth',4,'color',[0.7 0.7 0.7]);
        hold on; plot(XX,[0 0],'--','LineWidth',1,'color',[0.2 0.2 0.2]);
        hold on; plot([0 0],YY,'--','LineWidth',1,'color',[0.2 0.2 0.2]);
        if iplottype == 1 
            [~, pp_ttest] = ttest(XXX,YYY);
            [ pp_wil,  ~] = signrank(XXX,YYY);
            %
            XYdiag = [min([XX(1) YY(1)]),max([XX(2) YY(2)])];
            hold on; plot(XYdiag,XYdiag,'--','LineWidth',1,'color',[0.2 0.2 0.2]); 
        end
        plot(XXX, YYY, 'ko','MarkerSize',10,'LineWidth',1);
        Sigma_ell = cov(XXX, YYY);
        mu_ell(1) = mean(XXX);
        mu_ell(2) = mean(YYY);  
        % hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90); 
        if iplottype > 1 
            text(XX(1)+(abs((min(XX)))+abs((max(XX))))/20,YY(2)-(abs((min(YY)))+abs((max(YY))))/20,...
                {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
                 ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')]}, 'fontsize', 10);
        elseif iplottype == 1 
            text(XX(1)+(abs((min(XX)))+abs((max(XX))))/20,YY(2)-(abs((min(YY)))+abs((max(YY))))/20,...
                {['Spearman: r = ',num2str(RR_Spe,'%.2f'),', p = ',num2str(pp_Spe,'%1.1g')]; ...
                 ['Pearson: r = ',num2str(RR_Pea,'%.2f'),', p = ',num2str(pp_Pea,'%1.1g')];...
                 ['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                 ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                 ['mean \Delta = ',num2str(nanmean(YYY-XXX),'%.2f')];...
                 }, 'fontsize', 10);    
        end
        xlabel([xlabelname]);   
        ylabel([ylabelname]);  
        box off; axis([XX,YY]); axis square  
        set(gca,'fontsize',13) 
    end    

    % add monkey names
    axes('position',[.055 .70 .2 .05]); 
    hold on; h = text(0,0,{['monkey ',monkey_ana(1)]},'FontSize',14);axis off
    set(h,'Rotation',90);
end



end % for imonkey

% % % 
% % % 
% % % 
function plotErrorEllipse(mu_ell, Sigma_ell, p_ell)
    s = -2 * log(1 - p_ell);
    [V, D] = eig(Sigma_ell * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    plot(a(1, :) + mu_ell(1), a(2, :) + mu_ell(2), 'LineWidth',1.5, 'Color',[0.6 0.6 0.6] );
end


% % %
% % %
% % %
function plot_sigmoidcurve(cellname)
%
session = cellname(1:8);
readsession_TT;
%
filename = [dirroot,cellname,'_tuning'];
eval(['load ',filename])
table01_JC = tuning.JC.AB.table01;
table01_SO = abs(tuning.SO.ABA.table02);
%
filename = [dirroot,cellname,'_psyphycell'];
eval(['load ',filename])
% 
leg{1}='Task 2 (AB)';leg{2}='Task 2 (BA)'; leg{3}='Task 1';    
fitt = 'probit'; 
%
linew=3; dsize=70;
ptsord={'v';'^';'o'};
colord={[1. .2 .2];[ .2 .2 1.];[.5 .5 .5]};   
%
% SO
%separate non-forced choice trials and forced choices
table1mod_SO = table01_SO(logical(table01_SO(:,1) & table01_SO(:,2)),:);
forcedAtab_SO = table01_SO(logical(table01_SO(:,1) & ~table01_SO(:,2)),:);
forcedBtab_SO = table01_SO(logical(~table01_SO(:,1) & table01_SO(:,2)),:);
nfA = size(forcedAtab_SO,1);
nfB = size(forcedBtab_SO,1);
%
xx_SO = table1mod_SO(:,2)./table1mod_SO(:,1);
xx_SO = log(xx_SO);
yy_SO = table1mod_SO(:,3);
yy2_SO = table1mod_SO(:,3).*table1mod_SO(:,4);
%     
% seperate orders
xx3_SO = [];
for n=1:numel(yy_SO)
    if n<=(numel(yy_SO))/2; ordere='AB'; else ordere='BA'; end
    xx3_SO{n}=ordere;
end
xx3a_SO=strcmp(xx3_SO,'AB');
xx3b_SO=strcmp(xx3_SO,'BA');
xx3_SO=xx3a_SO-xx3b_SO;
%
warning off
Binosize_SO=table1mod_SO(:,4);
tbl=table(xx_SO,xx3_SO',yy2_SO,'VariableNames',{'offer','order','choice'});
[tbl idx]=sortrows(tbl,'order','descend');
Binosize_SO=Binosize_SO(idx);
mdl_SO=fitglm(tbl, 'choice ~ offer + order ','Distribution','binomial','link',fitt,'BinomialSize',Binosize_SO); % tell glmfit to use the binomial response
%
ord2=[1;-1]; x=[-2:0.025:3];
for ord=1:2
    for n=1:numel(x)
        order2(n)=ord2(ord);
    end
    tbl_pred{ord}=table(x',order2','VariableNames',{'offer','order'});

    [y_choicefit{ord}, ystd{ord}] = predict(mdl_SO,tbl_pred{ord});	
end
%
% JC
table1mod_JC = table01_JC(logical(table01_JC(:,1) & table01_JC(:,2)),:);
forcedAtab_JC = table01_JC(logical(table01_JC(:,1) & ~table01_JC(:,2)),:);
forcedBtab_JC = table01_JC(logical(~table01_JC(:,1) & table01_JC(:,2)),:);
nfA = size(forcedAtab_JC,1);
nfB = size(forcedBtab_JC,1);
%
xx_JC = table1mod_JC(:,2)./table1mod_JC(:,1);
xx_JC = log(xx_JC);
yy_JC = table1mod_JC(:,3);
yy2_JC = table1mod_JC(:,3).*table1mod_JC(:,4);
%
warning off
Binosize_JC=table1mod_JC(:,4);
tbl=table(xx_JC,yy2_JC,'VariableNames',{'offer','choice'});
mdl_JC=fitglm(tbl, 'choice ~ offer ','Distribution','binomial','link',fitt,'BinomialSize',Binosize_JC); % tell glmfit to use the binomial response
%
tbl_pred{3}=table(x','VariableNames',{'offer'});
[y_choicefit{3}, ystd{3}] = predict(mdl_JC,tbl_pred{3});		

%
% plot
lord=(numel(yy_SO)/2);    
%
hold on;
for ord=[3 1 2] % JC AB BA
    if ord < 3 % AB BA
        intval=0;
        hht(ord)=plot(x,y_choicefit{ord},'-','color',colord{ord},'linewidth',linew,'DisplayName',leg{ord});           
        sscat1 = scatter(xx_SO(lord*(ord-1)+1:lord*(ord)),yy_SO(lord*(ord-1)+1:lord*(ord)),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord},'DisplayName','off');
        sscat1.MarkerFaceAlpha = .4;
        %plot choice pattern + forced choices
        xx_forcedA = min(xx_SO)-log(1.3)*[1:nfA];	xx_forcedA = sort(xx_forcedA)';
        xx_forcedB = max(xx_SO)+log(1.3)*[1:nfB];	xx_forcedB = sort(xx_forcedB)';
        sscat2 = scatter(xx_forcedA,zeros(1,nfA),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord});
        sscat2.MarkerFaceAlpha = .4;
        sscat3 = scatter(xx_forcedB, ones(1,nfB),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord},'DisplayName','off');
        sscat3.MarkerFaceAlpha = .4;    
    else
        intval=0;
        hht(ord)=plot(x,y_choicefit{ord},'-','color',colord{ord},'linewidth',linew,'DisplayName',leg{ord});
        sscat1 = scatter(xx_JC,yy_JC,dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord});
        sscat1.MarkerFaceAlpha = .4;
        %plot choice pattern + forced choices
        xx_forcedA = min(xx_JC)-log(1.3)*[1:nfA];	xx_forcedA = sort(xx_forcedA)';
        xx_forcedB = max(xx_JC)+log(1.3)*[1:nfB];	xx_forcedB = sort(xx_forcedB)';
        sscat2 = scatter(xx_forcedA,zeros(1,nfA),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord});
        sscat2.MarkerFaceAlpha = .4;
        sscat3 = scatter(xx_forcedB, ones(1,nfB),dsize,ptsord{ord},'markerfacecolor',colord{ord},'markeredgecolor',colord{ord});
        sscat3.MarkerFaceAlpha = .4;                   
    end
end
%
legend 'show';
lgd = legend(hht([3,1,2]),'Location','northwest');
%
%cosmetics
set(gca,'xlim',[min(xx_SO)-1.1*log(1.3) max(xx_SO)+1.1*log(1.3)])
box off; axis square
[~,ind,~] = unique(xx_SO);	%remove doubles in xx
xxx = xx_SO(ind);   
% 
xlab = [];		%xlabels
for ifA = 1:nfA
    xlab{ifA} = [num2str(forcedAtab_SO(nfA-ifA+1,2)),':',num2str(forcedAtab_SO(nfA-ifA+1,1))];
end
for i = 1:size(xxx,1)
    xlab{nfA+i} = [num2str(table1mod_SO(ind(i),2)),':',num2str(table1mod_SO(ind(i),1))];
end
for ifB = 1:nfB
    xlab{nfA+i+ifB} = [num2str(forcedBtab_SO(nfB-ifB+1,2)),':',num2str(forcedBtab_SO(nfB-ifB+1,1))];
end
%
%add forced choices
xxx = [xx_forcedA;xxx;xx_forcedB];
set(gca,'fontsize',11)   
lgd.FontSize = 13;
set(gca,'xtick',xxx,'xticklabel',xlab)
set(gca,'ylim',[0,1]);
set(gca,'ytick',[0:.25:1],'yticklabel',[0:25:100]);
xtickangle(45);
xlabel('log(qB:qA)');
ylabel('B choice %');
title(cellname(1:8));
%
xpos=min(xx_SO)+(max(xx_SO)-min(xx_SO))*0.8;

text(xpos,.300,{['Task 1'];...
                ['\rho = ',	sprintf('%.2f',psyphycell.sigmoidfit.JC{3}(1))];...
                ['\eta = ',sprintf('%.2f',psyphycell.sigmoidfit.JC{2}(2))];[''];...
                ['Task 2'];...
                ['\rho = ',	sprintf('%.2f',psyphycell.sigmoidfit.SO{3}(1))];...
                ['\eta = ',sprintf('%.2f',psyphycell.sigmoidfit.SO{2}(2))];...
                ['\epsilon = ',sprintf('%.2f',(2*psyphycell.sigmoidfit.SO{3}(1)*psyphycell.sigmoidfit.SO{2}(3)/psyphycell.sigmoidfit.SO{2}(2)))]},'fontsize',11)
end






