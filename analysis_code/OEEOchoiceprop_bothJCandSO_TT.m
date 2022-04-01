% OEEOchoiceprop_bothJCandSO_TT.m
%
% This script examines the choice probability in/after offer2 time window in Task 2 (SO) in TT.
% This choice probability may reflect the working memory of decision and related to order bias or steepness or preference bias
% This script only examine chosen juice cells
% choice probability are calculated:
% 1. with all trials
% 2. without force choice
% 3. hard trials
% 4. easy trials

% author:   Feb 2021 - WS

close all
clearvars

brainarea = 'OFC'; % 'DLPFC', 'VLPFC' 'OFC'
monkey_ana = 'both'; % 'Gervinho', 'Juan', 'both'

dobhvlogit = ''; % '': probit; '_logit'; '_neworderbias'

mintrialnum = 100;

differentRho = '_differRho';    % '_differRho' or ''; do analysis based on the same rho in JC and SO or not 
twoTWinJCandSO = '';            % '' or ''; do fewer TW: 2 for JC and 2 for SO
doJCseq = '';                   % '_JCseq' or ''; do sequential JC parameters or not
nntrials = '_2nntrials';        % '_2nntrials' or '_3nntrials'
anovasetup = '_01pvallessTW_95nonzero'; % '_05pval_95nonzero'; '_01pvallessTW_95nonzero';
filename = ['TTcellist',doJCseq,'_',brainarea,'_both',differentRho,nntrials,anovasetup];
load(filename);

% cell analysis types - identify cells based on different conditions
JCSOclassifySep = 0; % if 1, load neurons that is classified seperately by JC or SO trials: 0-JCSOcellist
                        % only use 'SO' in this script

JCSOs = {'SO'};
nJCSO = size(JCSOs,2);       

alignments = {'off','cho'}; % 'off', 'cho', 'cue'
naligns = size(alignments,2);

cellclassnames = {'CJ'};
classnamenums = [4];
nclassnames = size(cellclassnames,2);

%
if JCSOclassifySep
    JCclasses = [JCcellist.infocell.subclass]';
    JCslopes  = [JCcellist.infocell.slopesign]';
    SOclasses = [SOcellist.infocell.subclass]';
    SOslopes  = [SOcellist.infocell.slopesign]';
    cellnames = {JCcellist.infocell.cellname}';
elseif ~JCSOclassifySep
    JCclasses = [JCSOcellist.infocell.subclass]';
    JCslopes  = [JCSOcellist.infocell.slopesign]';
    SOclasses = [JCSOcellist.infocell.subclass]';
    SOslopes  = [JCSOcellist.infocell.slopesign]';
    cellnames = {JCSOcellist.infocell.cellname}';
end


% % % % % % 
% remove outlier based on steepness across all sessions (keep consistent with behavioral measurement)
removesessions = 0; % remove based on dynamic range and saturation
dosteepout = 0; 
% Gervinho
if removesessions
    filename = ['pop_behav_ana_summary_Gervinho_removesessions',dobhvlogit];
elseif ~removesessions
    filename = ['pop_behav_ana_summary_Gervinho',dobhvlogit];
end
eval(['load ', filename])
quantiles_JC = quantile(steepness_all(:,1),[0.25 0.5 0.75]);
IQR_JC = quantiles_JC(3) - quantiles_JC(1);
steepoutlier_JC_G = [quantiles_JC(1)-1.5*IQR_JC, quantiles_JC(3)+1.5*IQR_JC];
quantiles_SO = quantile(steepness_all(:,2),[0.25 0.5 0.75]);
IQR_SO = quantiles_SO(3) - quantiles_SO(1);
steepoutlier_SO_G = [quantiles_SO(1)-1.5*IQR_SO, quantiles_SO(3)+1.5*IQR_SO];
% Juan
if removesessions
    filename = ['pop_behav_ana_summary_Juan_removesessions',dobhvlogit];
elseif ~removesessions
    filename = ['pop_behav_ana_summary_Juan',dobhvlogit];
end
eval(['load ', filename])
quantiles_JC = quantile(steepness_all(:,1),[0.25 0.5 0.75]);
IQR_JC = quantiles_JC(3) - quantiles_JC(1);
steepoutlier_JC_J = [quantiles_JC(1)-1.5*IQR_JC, quantiles_JC(3)+1.5*IQR_JC];
quantiles_SO = quantile(steepness_all(:,2),[0.25 0.5 0.75]);
IQR_SO = quantiles_SO(3) - quantiles_SO(1);
steepoutlier_SO_J = [quantiles_SO(1)-1.5*IQR_SO, quantiles_SO(3)+1.5*IQR_SO];
% % % % % % 


%
for iclass = 1:nclassnames % nclassnames == 1; CJ only
    classnamenum = classnamenums(iclass);
    ind_JCiclass = JCclasses==classnamenum;
    ind_SOiclass = SOclasses==classnamenum;
    ind_iclass = ind_JCiclass & ind_SOiclass; % |: pass the classification in either JC or SO; &: pass the classification in both JC and SO
    % ind_iclass = ind_SOiclass;
    %
    ncells = sum(ind_iclass);
    cellnames_iclass = cellnames(ind_iclass);
    cellnum_iclass   = find(ind_iclass==1);
    JCslopes_iclass  = JCslopes(ind_iclass); % 1: CJA, -1 CJB
    SOslopes_iclass  = SOslopes(ind_iclass); % 1: CJA, -1 CJB

    
    for iJCSO = 1:nJCSO
        JCorSO = JCSOs{iJCSO};
        [~, ~,  timewindows_off, timewindows_cho, ~, timewindows_cue, ~] = get_timewindows_TT(JCorSO);
    
        for ialign = 1:naligns
            alignment = alignments{ialign};         
            eval(['ntwins = size(timewindows_',alignment,',2);'])
            ncells_ana = 0;
            % %
            timepoints = [];
            monkeyname_list = [];
            cellnames_list = [];
            tuningslopes_list = [];
            % 
            rho_all = [];
            steepness_all = [];    
            orderbias_all = [];
            %
            corr_VOandVE = [];
            % % 
            CP_in_OE_movTW = [];
            CP_in_EO_movTW = []; 
            CPpval_in_OE_movTW = [];
            CPpval_in_EO_movTW = [];
            % %
            FR_in_OE_movTW = [];
            FR_in_EO_movTW = [];
            %
            FR_in_OE_movTW_noforced = [];
            FR_in_EO_movTW_noforced = [];
            %
            FR_in_OE_movTW_hard = [];
            FR_in_EO_movTW_hard = [];
            %
            FR_in_OE_movTW_easy = [];
            FR_in_EO_movTW_easy = [];
            % %
            CP_in_OE_criTW = [];
            CP_in_EO_criTW = [];
            CPpval_in_OE_criTW = [];
            CPpval_in_EO_criTW = [];
            %
            CP_in_OE_criTW_noforced = [];
            CP_in_EO_criTW_noforced = [];
            CPpval_in_OE_criTW_noforced = [];
            CPpval_in_EO_criTW_noforced = [];
            %
            CP_in_OE_criTW_hard = [];
            CP_in_EO_criTW_hard = [];
            CPpval_in_OE_criTW_hard = [];
            CPpval_in_EO_criTW_hard = [];
            %
            CP_in_OE_criTW_easy = [];
            CP_in_EO_criTW_easy = [];
            CPpval_in_OE_criTW_easy = [];
            CPpval_in_EO_criTW_easy = [];
            % %
            CP_in_OEEO_criTW = [];
            CPpval_in_OEEO_criTW = [];
            %
            CP_in_OEEO_criTW_noforced = [];
            CPpval_in_OEEO_criTW_noforced = [];
            %
            CP_in_OEEO_criTW_hard = [];
            CPpval_in_OEEO_criTW_hard = [];
            %
            CP_in_OEEO_criTW_easy = [];
            CPpval_in_OEEO_criTW_easy = [];
            
            
            try
                % exnovo
                filename = ['OEEOchoiceprop_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,'',dobhvlogit];                   
                load(filename); 
            
            catch
            
                for icell = 1:ncells
                    cellname = cellnames_iclass{icell};
                    session = cellname(1:8);
                    readsession_TT;

                    %%%%%%%%%%%
                    if 1
                    % remove sessions accordingly 
                    % skip sessions with too few trials
                    filename = [dirroot,cellname,'_data'];
                    eval(['load ',filename])
                    clear celldata celldataerror trace trialRecord
                    if size(goodTrials_JC,1)<mintrialnum | size(goodTrials_SO,1)<mintrialnum
                        continue
                    end   
                    if psyphydata(find(psyphydata(:,3)==33,1),1) - psyphydata(find(psyphydata(:,3)==31,1),1) > 1150
                        continue
                    end
                    end
                    %%%%%%%%%%%
                    
                    ncells_ana = ncells_ana+1;
                    %
                    monkeyname_list{ncells_ana} = cellname(1);
                    cellnames_list{ncells_ana} = cellname;
            
                    %
                    JCslope_icell = JCslopes_iclass(icell);
                    SOslope_icell = SOslopes_iclass(icell);
                    tuningslopes_list(ncells_ana) = SOslope_icell;
                    
                    %
                    filename = [dirroot,cellname,'_data'];
                    eval(['load ',filename])
                    clear celldataerror trace trialRecord
                    ind = ismember(psyphydata(:,3),[43,44,45]);
                    psyphydata(ind,3) = 50;
                    
                    
                    %
                    try
                    % exnovo
                    filename = [dirroot, cellname,'_tuning'];
                    eval(['load ',filename])
                    catch
                    tuning_TT(cellname);    % CHANGE IT ALIGNED TO OFFER 2 IN SO
                    filename = [dirroot, cellname,'_tuning'];
                    eval(['load ',filename])
                    end
                    
                   
                    % 
                    % % load behavior
                    try 
                        dummy
                        filename = [dirroot,cellname,'_psyphycell'];
                        eval(['load ',filename])                                    
                        if isempty(differentRho)
                            rho_JC = (psyphycell.sigmoidfit.JC{3}(1)+psyphycell.sigmoidfit.SO{3}(1))/2;
                            rho_SO = (psyphycell.sigmoidfit.JC{3}(1)+psyphycell.sigmoidfit.SO{3}(1))/2;
                        elseif ~isempty(differentRho)
                            rho_JC = psyphycell.sigmoidfit.JC{3}(1);                        
                            rho_SO = psyphycell.sigmoidfit.SO{3}(1); 
                        end
                        steepness_JC = psyphycell.sigmoidfit.JC{2}(2);
                        steepness_SO = psyphycell.sigmoidfit.SO{2}(2);
                        orderbias_SO = -2.*rhos_inSO(icell,:).*psyphycell.sigmoidfit.SO{2}(3)./psyphycell.sigmoidfit.SO{2}(2);    
                    catch
                        warning off
                        % [psyphycell] = sigmoidfit_TT_OrdChHyst([cellname],'probit','Only',1); % 'Only': SO on SO or JC on JC; 'Both': SO and JC on either SO or JC
                        if isequal(dobhvlogit,'_logit')
                            [psyphycell] = sigmoidfit_TT_OrdChHyst([cellname],'logit','Both','log',1);
                        else
                            [psyphycell] = sigmoidfit_TT_OrdChHyst([cellname],'probit','Both','log',1);
                        end
                        relvalue_JC = exp(-(psyphycell.JC.NonChHyst.sigmoidfit.beta(1))/(psyphycell.JC.NonChHyst.sigmoidfit.beta(2)));
                        relvalue_SO = exp(-(psyphycell.SO.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SO.NonChHyst.sigmoidfit.beta(2)));                                    
                        if isempty(differentRho)
                            rho_JC = (relvalue_SO + relvalue_JC)/2;
                            rho_SO = (relvalue_SO + relvalue_JC)/2;
                        elseif ~isempty(differentRho)
                            rho_JC = relvalue_JC;
                            rho_SO = relvalue_SO;
                        end 
                        steepness_JC = psyphycell.JC.NonChHyst.sigmoidfit.beta(2);
                        steepness_SO = psyphycell.SO.NonChHyst.sigmoidfit.beta(2);     
                        if isequal(dobhvlogit, '_neworderbias')
                            rho_inAB = exp(-(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(2)));
                            rho_inBA = exp(-(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(2)));
                            orderbias_inSO(icell,:) = rho_inBA - rho_inAB;
                        else
                            orderbias_SO = -2.*rho_SO.*psyphycell.SO.NonChHyst.sigmoidfit.beta(3)./psyphycell.SO.NonChHyst.sigmoidfit.beta(2);
                        % orderbias_SO = -2.*rho_SO.*psyphycell.SO.OrdChHyst.sigmoidfit.beta(3)./psyphycell.SO.OrdChHyst.sigmoidfit.beta(2);
                        end                       
                    end
                    %
                    steepness_all.JC(:,ncells_ana) = steepness_JC;
                    rho_all.JC(:,ncells_ana) = rho_JC;
                    %
                    steepness_all.SO(:,ncells_ana) = steepness_SO;
                    rho_all.SO(:,ncells_ana) = rho_SO;
                    orderbias_all.SO(:,ncells_ana) = orderbias_SO;
                    
                    
                    %
                    % load neuron data 
                    % JC trials
                    if isequal(JCorSO,'JC')
                        pairname = 'AB';
                        eval(['superneuract=tuning.',JCorSO,'.',pairname,'.superneuract.bytrial;']) 
                       
                    % SO trials
                    elseif isequal(JCorSO,'SO')
                        pairname = 'ABA';
                        eval(['superneuract=tuning.',JCorSO,'.',pairname,'.superneuract.bytrial;'])
                    end % if JC or SO    
                  
                    
                    %     
                    % moving time window; choice probability 
                    allneuroFR_icell = [];
                    for itwin = 1:ntwins
                        eval(['timepoints(itwin,1) = mean(timewindows_',alignment,'{itwin}{3});'])
                        %
                        eval(['neuract = superneuract.twds_',alignment,num2str(itwin),';'])     
                        % 
                        % % remove force choice
                        % ind_forced = neuract(:,2)==0 | neuract(:,3)==0;
                        % neuract(ind_forced,:) = [];
                        
                        
                        % JC trials
                        if isequal(JCorSO,'JC')
                            % under construction
                        %       
                        % SO trials   
                        elseif isequal(JCorSO,'SO')
                            %
                            ind_forced = neuract(:,2)==0 | neuract(:,3)==0;
                            %
                            % hard and easy
                            % hard trials (choice of B is 10%-90%) and easy trials
                            % define hard and easy trials
                            if isequal(monkeyname_list{ncells_ana},'G')
                                eps = 0.1;
                            elseif isequal(monkeyname_list{ncells_ana},'J')
                                eps = 0.1; % 0.2
                            else
                                eps = 0.1;
                            end
                            atleast_nntrials = 2;
                            %
                            OVA_SO = neuract(:,2);
                            OVB_SO = neuract(:,3);                           
                            CJ_SO = neuract(:,5);  % -1: B
                            ord_SO = neuract(:,6); % -1: BA                               
                            % average over trial type: [OVA, OVB, ord]
                            % add other parameters into trial type table
                            % [OVA, OVB, ord, #trials, propB]
                            newTrialTypes_all = [OVA_SO, OVB_SO, ord_SO];
                            uni_newTriTypes = unique(newTrialTypes_all,'rows');
                            nnewTriTypes = size(uni_newTriTypes,1);
                            TrialTypes_tbl = uni_newTriTypes;
                            for inewTriType = 1:nnewTriTypes
                                newTriType = uni_newTriTypes(inewTriType,:);
                                ind_newTriType = ismember(newTrialTypes_all, newTriType,'rows');
                                %
                                trialnum = sum(ind_newTriType);
                                TrialTypes_tbl(inewTriType,4) = trialnum;
                                %
                                trialnum_chB = sum(CJ_SO(ind_newTriType)==-1);
                                TrialTypes_tbl(inewTriType,5) = trialnum_chB./trialnum;                                
                            end
                            % 
                            ind_hardtypes = TrialTypes_tbl(:,5)<1-eps & TrialTypes_tbl(:,5)>eps;
                            TrialTypes_tbl_hard = TrialTypes_tbl(ind_hardtypes,:);
                            %
                            ind_hard = ismember(neuract(:,[2,3,6]),TrialTypes_tbl_hard(:,1:3),'rows');
                            ind_easy = ~ind_hard & ~ind_forced; 
                            
                            %
                            % OE hard trials
                            orderABBA = neuract(:,6); % -1: BA 1: AB
                            neurFR = neuract(:,7);
                            if SOslope_icell == 1 % CJA
                                ind_BA = orderABBA == -1 & ind_hard;
                                neuract_OE = neuract(ind_BA,:);
                                Echosen = neuract_OE(:,5) == 1;                                
                            elseif SOslope_icell == -1 %CJB
                                ind_AB = orderABBA == 1 & ind_hard;
                                neuract_OE = neuract(ind_AB,:);
                                Echosen = neuract_OE(:,5) == -1; 
                            end 
                            trialset = {};
                            trialset(Echosen==1) = {'Ehard'};
                            trialset(Echosen==0) = {'Ohard'};
                            fr_act = neuract_OE(:,7);
                            try
                                [~,~,~,auc] = perfcurve(trialset,fr_act,'Ehard','nboot',100,'alpha',.05);
                                %
                                CP_in_OE_movTW(itwin,ncells_ana) = auc(1);
                                CPpval_in_OE_movTW(itwin,ncells_ana) = prod(auc(2:3)-.5)>0;
                            catch
                                CP_in_OE_movTW(itwin,ncells_ana) = nan;
                                CPpval_in_OE_movTW(itwin,ncells_ana) = 0;
                            end
                            %
                            %
                            % EO hard trials
                            orderABBA = neuract(:,6); % -1: BA 1: AB
                            if SOslope_icell == 1 % CJA
                                ind_AB = orderABBA == 1 & ind_hard;
                                neuract_EO = neuract(ind_AB,:);
                                Echosen = neuract_EO(:,5) == 1;                                
                            elseif SOslope_icell == -1 %CJB
                                ind_BA = orderABBA == -1 & ind_hard;
                                neuract_EO = neuract(ind_BA,:);
                                Echosen = neuract_EO(:,5) == -1; 
                            end 
                            trialset = {};
                            trialset(Echosen==1) = {'Ehard'};
                            trialset(Echosen==0) = {'Ohard'};
                            fr_act = neuract_EO(:,7);
                            try
                                [~,~,~,auc] = perfcurve(trialset,fr_act,'Ehard','nboot',100,'alpha',.05);
                                %
                                CP_in_EO_movTW(itwin,ncells_ana) = auc(1);
                                CPpval_in_EO_movTW(itwin,ncells_ana) = prod(auc(2:3)-.5)>0;
                            catch
                                CP_in_EO_movTW(itwin,ncells_ana) = nan;
                                CPpval_in_EO_movTW(itwin,ncells_ana) = 0;
                            end
                        end %if else JC SO  
                        allneuroFR_icell(:,itwin) = neurFR; % each row: each trial; each colomn: each time window 
                    end %for itwin        

                    
                    % %
                    % average firing rate across critical time window                    
                    % 
                    if isequal(JCorSO,'JC')
                        %
                        % UNDER CONSTRUCTION!!
                    elseif isequal(JCorSO,'SO')
                        if isequal(alignment,'off')     
                            alignpoint = [30 30];
                            % TIME WINDOWS: 0-500ms offer1; 500-1000ms interoffer; 1000-1500ms offer2; 1500ms-2000ms wait                        
                            critTWs = [-250,    0;   % before offer1 on as control
                                       1250, 1500;   % later during offer2 onset
                                       1500, 1750;   % early during waiting onset 
                                       1750, 2000;   % later during waiting onset
                                       2000, 2250;   % after waiting
                                       ]; 
                            critTWnames = {'beforeofferon','lateoffer2','earlywait','latewait','afterwait'};
                            ncriTWs = length(critTWnames);
                            %
                        elseif isequal(alignment,'cho')       
                            alignpoint = [50 50];
                            % TIME WINDOWS: 0ms: choice                 
                            critTWs = [-500  -250;   % early pre-juice
                                       -250,    0;   % pre-juice
                                          0,  250;   % post-juice
                                        250,  500;   % late post-juice
                                       ]; 
                            critTWnames = {'earlyprejuice', 'prejuice', 'postjuice','latepostjuice'};
                            ncriTWs = length(critTWnames);
                            %
                        end % isequal 'off', 'cho'
                        %
                        for icriTW = 1:ncriTWs
                            %
                            critTW = critTWs(icriTW,:);
                            critTWname = critTWnames{icriTW};
                            timewindow_new = {critTWname,alignpoint,critTW};
                            %    
                            ind_TW = timepoints>=critTW(1) & timepoints<=critTW(2);
                            averCritFR = nanmean(allneuroFR_icell(:,ind_TW),2); % each row: each trial; one colomn: average across critical TW
                            % 
                            % 
                            orderABBA = neuract(:,6); % -1: BA 1: AB
                            ind_forced = neuract(:,2)==0 | neuract(:,3)==0;
                            
                            % 
                            % all trials
                            % OE trials all
                            if SOslope_icell == 1 % CJA
                                ind_BA = orderABBA == -1;
                                ind_chosenE = ind_BA & neuract(:,5) ==  1;  
                                ind_chosenO = ind_BA & neuract(:,5) == -1; 
                                %
                                neuract_OE = neuract(ind_BA,:);
                                Echosen = neuract_OE(:,5) == 1;    
                                averCritFR_OE = averCritFR(ind_BA,:);
                                FR_in_OE = allneuroFR_icell(ind_BA,:);
                            elseif SOslope_icell == -1 %CJB
                                ind_AB = orderABBA == 1;
                                ind_chosenE = ind_AB & neuract(:,5) == -1;  
                                ind_chosenO = ind_AB & neuract(:,5) ==  1; 
                                %
                                neuract_OE = neuract(ind_AB,:);
                                Echosen = neuract_OE(:,5) == -1; 
                                averCritFR_OE = averCritFR(ind_AB,:);
                                FR_in_OE = allneuroFR_icell(ind_AB,:);
                            end 
                            hits_chosenE = neuract(ind_chosenE,1);
                            hits_chosenO = neuract(ind_chosenO,1);
                            neuract_chosenE = actmake_TT(celldata, psyphydata, hits_chosenE, timewindow_new);
                            neuract_chosenO = actmake_TT(celldata, psyphydata, hits_chosenO, timewindow_new);                              
                            Echosen_ROC = [ones(size(neuract_chosenE)); zeros(size(neuract_chosenO))];
                            fr_act = [neuract_chosenE; neuract_chosenO];
                            %
                            trialset = {};
                            trialset(Echosen_ROC == 1) = {'Eall'};
                            trialset(Echosen_ROC == 0) = {'Oall'};
                            [~,~,~,auc] = perfcurve(trialset,fr_act,'Eall','nboot',100,'alpha',.05);
                            %
                            eval(['CP_in_OE_criTW.',critTWname,'(:,ncells_ana) = auc(1);'])
                            eval(['CPpval_in_OE_criTW.',critTWname,'(:,ncells_ana) = prod(auc(2:3)-.5)>0;'])
                            %                           
                            FR_in_OE_movTW.Echosen(:,ncells_ana) = nanmean(FR_in_OE(Echosen==1,:));
                            FR_in_OE_movTW.Ochosen(:,ncells_ana) = nanmean(FR_in_OE(Echosen==0,:));
                            %
                            % EO trials all
                            if SOslope_icell == 1 % CJA
                                ind_AB = orderABBA == 1;
                                ind_chosenE = ind_AB & neuract(:,5) ==  1;  
                                ind_chosenO = ind_AB & neuract(:,5) == -1;  
                                %
                                neuract_EO = neuract(ind_AB,:);
                                Echosen = neuract_EO(:,5) == 1;    
                                averCritFR_EO = averCritFR(ind_AB,:);
                                FR_in_EO = allneuroFR_icell(ind_AB,:);
                            elseif SOslope_icell == -1 %CJB
                                ind_BA = orderABBA == -1;
                                ind_chosenE = ind_BA & neuract(:,5) == -1;  
                                ind_chosenO = ind_BA & neuract(:,5) ==  1; 
                                %
                                neuract_EO = neuract(ind_BA,:);
                                Echosen = neuract_EO(:,5) == -1; 
                                averCritFR_EO = averCritFR(ind_BA,:);
                                FR_in_EO = allneuroFR_icell(ind_BA,:);
                            end 
                            hits_chosenE = neuract(ind_chosenE,1);
                            hits_chosenO = neuract(ind_chosenO,1);
                            neuract_chosenE = actmake_TT(celldata, psyphydata, hits_chosenE, timewindow_new);
                            neuract_chosenO = actmake_TT(celldata, psyphydata, hits_chosenO, timewindow_new);                              
                            Echosen_ROC = [ones(size(neuract_chosenE)); zeros(size(neuract_chosenO))];
                            fr_act = [neuract_chosenE; neuract_chosenO];
                            %
                            trialset = {};
                            trialset(Echosen_ROC == 1) = {'Eall'};
                            trialset(Echosen_ROC == 0) = {'Oall'};
                            [~,~,~,auc] = perfcurve(trialset,fr_act,'Eall','nboot',100,'alpha',.05);
                            %
                            eval(['CP_in_EO_criTW.',critTWname,'(:,ncells_ana) = auc(1);'])
                            eval(['CPpval_in_EO_criTW.',critTWname,'(:,ncells_ana) = prod(auc(2:3)-.5)>0;'])
                            %                           
                            FR_in_EO_movTW.Echosen(:,ncells_ana) = nanmean(FR_in_EO(Echosen==1,:));
                            FR_in_EO_movTW.Ochosen(:,ncells_ana) = nanmean(FR_in_EO(Echosen==0,:));
                                                       
                                                      
                            % 
                            % all trials without forced choice
                            % OE trials without forced choice
                            if SOslope_icell == 1 % CJA
                                ind_BA = orderABBA ==-1 & ~ind_forced;
                                ind_chosenE = ind_BA & neuract(:,5) ==  1;  
                                ind_chosenO = ind_BA & neuract(:,5) == -1; 
                                %
                                neuract_OE = neuract(ind_BA,:);
                                Echosen = neuract_OE(:,5) == 1;    
                                averCritFR_OE = averCritFR(ind_BA,:);
                                FR_in_OE = allneuroFR_icell(ind_BA,:);
                            elseif SOslope_icell == -1  %CJB
                                ind_AB = orderABBA == 1 & ~ind_forced;
                                ind_chosenE = ind_AB & neuract(:,5) == -1;  
                                ind_chosenO = ind_AB & neuract(:,5) ==  1; 
                                %
                                neuract_OE = neuract(ind_AB,:);
                                Echosen = neuract_OE(:,5) == -1; 
                                averCritFR_OE = averCritFR(ind_AB,:);
                                FR_in_OE = allneuroFR_icell(ind_AB,:);
                            end 
                            hits_chosenE = neuract(ind_chosenE,1);
                            hits_chosenO = neuract(ind_chosenO,1);
                            neuract_chosenE = actmake_TT(celldata, psyphydata, hits_chosenE, timewindow_new);
                            neuract_chosenO = actmake_TT(celldata, psyphydata, hits_chosenO, timewindow_new);                              
                            Echosen_ROC = [ones(size(neuract_chosenE)); zeros(size(neuract_chosenO))];
                            fr_act = [neuract_chosenE; neuract_chosenO];
                            %
                            trialset = {};
                            trialset(Echosen_ROC==1) = {'Eallnoforced'};
                            trialset(Echosen_ROC==0) = {'Oallnoforced'};
                            [~,~,~,auc] = perfcurve(trialset,fr_act,'Eallnoforced','nboot',100,'alpha',.05);
                            %
                            eval(['CP_in_OE_criTW_noforced.',critTWname,'(:,ncells_ana) = auc(1);'])
                            eval(['CPpval_in_OE_criTW_noforced.',critTWname,'(:,ncells_ana) = prod(auc(2:3)-.5)>0;'])
                            %                           
                            FR_in_OE_movTW_noforced.Echosen(:,ncells_ana) = nanmean(FR_in_OE(Echosen==1,:));
                            FR_in_OE_movTW_noforced.Ochosen(:,ncells_ana) = nanmean(FR_in_OE(Echosen==0,:));
                            %
                            % EO trials without forced choice
                            if SOslope_icell == 1 % CJA
                                ind_AB = orderABBA == 1 & ~ind_forced;
                                ind_chosenE = ind_AB & neuract(:,5) ==  1;  
                                ind_chosenO = ind_AB & neuract(:,5) == -1;  
                                %
                                neuract_EO = neuract(ind_AB,:);
                                Echosen = neuract_EO(:,5) == 1;    
                                averCritFR_EO = averCritFR(ind_AB,:);
                                FR_in_EO = allneuroFR_icell(ind_AB,:);
                            elseif SOslope_icell == -1 %CJB
                                ind_BA = orderABBA ==-1  & ~ind_forced;
                                ind_chosenE = ind_BA & neuract(:,5) == -1;  
                                ind_chosenO = ind_BA & neuract(:,5) ==  1; 
                                %
                                neuract_EO = neuract(ind_BA,:);
                                Echosen = neuract_EO(:,5) == -1; 
                                averCritFR_EO = averCritFR(ind_BA,:);
                                FR_in_EO = allneuroFR_icell(ind_BA,:);
                            end 
                            hits_chosenE = neuract(ind_chosenE,1);
                            hits_chosenO = neuract(ind_chosenO,1);
                            neuract_chosenE = actmake_TT(celldata, psyphydata, hits_chosenE, timewindow_new);
                            neuract_chosenO = actmake_TT(celldata, psyphydata, hits_chosenO, timewindow_new);                              
                            Echosen_ROC = [ones(size(neuract_chosenE)); zeros(size(neuract_chosenO))];
                            fr_act = [neuract_chosenE; neuract_chosenO];
                            %
                            trialset = {};
                            trialset(Echosen_ROC==1) = {'Eallnoforced'};
                            trialset(Echosen_ROC==0) = {'Oallnoforced'};
                            [~,~,~,auc] = perfcurve(trialset,fr_act,'Eallnoforced','nboot',100,'alpha',.05);
                            %
                            eval(['CP_in_EO_criTW_noforced.',critTWname,'(:,ncells_ana) = auc(1);'])
                            eval(['CPpval_in_EO_criTW_noforced.',critTWname,'(:,ncells_ana) = prod(auc(2:3)-.5)>0;'])
                            %
                            FR_in_EO_movTW_noforced.Echosen(:,ncells_ana) = nanmean(FR_in_EO(Echosen==1,:));
                            FR_in_EO_movTW_noforced.Ochosen(:,ncells_ana) = nanmean(FR_in_EO(Echosen==0,:));
                            
                            
                            % 
                            % hard trials (choice of B is 10%-90%) and easy trials
                            % define hard and easy trials
                            if isequal(monkeyname_list{ncells_ana},'G')
                                eps = 0.1;
                            elseif isequal(monkeyname_list{ncells_ana},'J')
                                eps = 0.1; % 0.2
                            else
                                eps = 0.1;
                            end
                            atleast_nntrials = 2;
                            %
                            OVA_SO = neuract(:,2);
                            OVB_SO = neuract(:,3);                           
                            CJ_SO = neuract(:,5);  % -1: B
                            ord_SO = neuract(:,6); % -1: BA                               
                            % average over trial type: [OVA, OVB, ord]
                            % add other parameters into trial type table
                            % [OVA, OVB, ord, #trials, propB]
                            newTrialTypes_all = [OVA_SO, OVB_SO, ord_SO];
                            uni_newTriTypes = unique(newTrialTypes_all,'rows');
                            nnewTriTypes = size(uni_newTriTypes,1);
                            TrialTypes_tbl = uni_newTriTypes;
                            for inewTriType = 1:nnewTriTypes
                                newTriType = uni_newTriTypes(inewTriType,:);
                                ind_newTriType = ismember(newTrialTypes_all, newTriType,'rows');
                                %
                                trialnum = sum(ind_newTriType);
                                TrialTypes_tbl(inewTriType,4) = trialnum;
                                %
                                trialnum_chB = sum(CJ_SO(ind_newTriType)==-1);
                                TrialTypes_tbl(inewTriType,5) = trialnum_chB./trialnum;                                
                            end
                            % 
                            ind_hardtypes = TrialTypes_tbl(:,5)<1-eps & TrialTypes_tbl(:,5)>eps;
                            TrialTypes_tbl_hard = TrialTypes_tbl(ind_hardtypes,:);
                            %
                            ind_hard = ismember(neuract(:,[2,3,6]),TrialTypes_tbl_hard(:,1:3),'rows');
                            ind_easy = ~ind_hard & ~ind_forced;    
                            
                            %
                            %
                            try 
                                % OE trials hard
                                if SOslope_icell == 1 % CJA
                                    ind_BA = orderABBA ==-1 & ind_hard;
                                    ind_chosenE = ind_BA & neuract(:,5) ==  1;  
                                    ind_chosenO = ind_BA & neuract(:,5) == -1; 
                                    %
                                    neuract_OE = neuract(ind_BA,:);
                                    Echosen = neuract_OE(:,5) == 1;    
                                    averCritFR_OE = averCritFR(ind_BA,:);
                                    FR_in_OE = allneuroFR_icell(ind_BA,:);
                                elseif SOslope_icell == -1  %CJB
                                    ind_AB = orderABBA == 1 & ind_hard;
                                    ind_chosenE = ind_AB & neuract(:,5) == -1;  
                                    ind_chosenO = ind_AB & neuract(:,5) ==  1; 
                                    %
                                    neuract_OE = neuract(ind_AB,:);
                                    Echosen = neuract_OE(:,5) == -1; 
                                    averCritFR_OE = averCritFR(ind_AB,:);
                                    FR_in_OE = allneuroFR_icell(ind_AB,:);
                                end 
                                hits_chosenE = neuract(ind_chosenE,1);
                                hits_chosenO = neuract(ind_chosenO,1);
                                neuract_chosenE = actmake_TT(celldata, psyphydata, hits_chosenE, timewindow_new);
                                neuract_chosenO = actmake_TT(celldata, psyphydata, hits_chosenO, timewindow_new);                              
                                Echosen_ROC = [ones(size(neuract_chosenE)); zeros(size(neuract_chosenO))];
                                fr_act = [neuract_chosenE; neuract_chosenO];
                                %
                                trialset = {};
                                trialset(Echosen_ROC==1) = {'Ehard'};
                                trialset(Echosen_ROC==0) = {'Ohard'};
                                [~,~,~,auc] = perfcurve(trialset,fr_act,'Ehard','nboot',100,'alpha',.05);
                                %
                                eval(['CP_in_OE_criTW_hard.',critTWname,'(:,ncells_ana) = auc(1);'])
                                eval(['CPpval_in_OE_criTW_hard.',critTWname,'(:,ncells_ana) = prod(auc(2:3)-.5)>0;'])
                                %                           
                                FR_in_OE_movTW_hard.Echosen(:,ncells_ana) = nanmean(FR_in_OE(Echosen==1,:));
                                FR_in_OE_movTW_hard.Ochosen(:,ncells_ana) = nanmean(FR_in_OE(Echosen==0,:));
                                %
                                % EO trials hard
                                if SOslope_icell == 1 % CJA
                                    ind_AB = orderABBA == 1 & ind_hard;
                                    ind_chosenE = ind_AB & neuract(:,5) ==  1;  
                                    ind_chosenO = ind_AB & neuract(:,5) == -1;  
                                    %
                                    neuract_EO = neuract(ind_AB,:);
                                    Echosen = neuract_EO(:,5) == 1;    
                                    averCritFR_EO = averCritFR(ind_AB,:);
                                    FR_in_EO = allneuroFR_icell(ind_AB,:);
                                elseif SOslope_icell == -1 %CJB
                                    ind_BA = orderABBA ==-1  & ind_hard;
                                    ind_chosenE = ind_BA & neuract(:,5) == -1;  
                                    ind_chosenO = ind_BA & neuract(:,5) ==  1; 
                                    %
                                    neuract_EO = neuract(ind_BA,:);
                                    Echosen = neuract_EO(:,5) == -1; 
                                    averCritFR_EO = averCritFR(ind_BA,:);
                                    FR_in_EO = allneuroFR_icell(ind_BA,:);
                                end 
                                hits_chosenE = neuract(ind_chosenE,1);
                                hits_chosenO = neuract(ind_chosenO,1);
                                neuract_chosenE = actmake_TT(celldata, psyphydata, hits_chosenE, timewindow_new);
                                neuract_chosenO = actmake_TT(celldata, psyphydata, hits_chosenO, timewindow_new);                              
                                Echosen_ROC = [ones(size(neuract_chosenE)); zeros(size(neuract_chosenO))];
                                fr_act = [neuract_chosenE; neuract_chosenO];
                                %
                                trialset = {};
                                trialset(Echosen_ROC==1) = {'Ehard'};
                                trialset(Echosen_ROC==0) = {'Ohard'};
                                [~,~,~,auc] = perfcurve(trialset,fr_act,'Ehard','nboot',100,'alpha',.05);
                                %
                                eval(['CP_in_EO_criTW_hard.',critTWname,'(:,ncells_ana) = auc(1);'])
                                eval(['CPpval_in_EO_criTW_hard.',critTWname,'(:,ncells_ana) = prod(auc(2:3)-.5)>0;'])
                                %
                                FR_in_EO_movTW_hard.Echosen(:,ncells_ana) = nanmean(FR_in_EO(Echosen==1,:));
                                FR_in_EO_movTW_hard.Ochosen(:,ncells_ana) = nanmean(FR_in_EO(Echosen==0,:));
                                
                                %
                                %
                                % OE trials easy
                                if SOslope_icell == 1 % CJA
                                    ind_BA = orderABBA ==-1 & ind_easy;
                                    ind_chosenE = ind_BA & neuract(:,5) ==  1;  
                                    ind_chosenO = ind_BA & neuract(:,5) == -1; 
                                    %
                                    neuract_OE = neuract(ind_BA,:);
                                    Echosen = neuract_OE(:,5) == 1;    
                                    averCritFR_OE = averCritFR(ind_BA,:);
                                    FR_in_OE = allneuroFR_icell(ind_BA,:);
                                elseif SOslope_icell == -1  %CJB
                                    ind_AB = orderABBA == 1 & ind_easy;
                                    ind_chosenE = ind_AB & neuract(:,5) == -1;  
                                    ind_chosenO = ind_AB & neuract(:,5) ==  1; 
                                    %
                                    neuract_OE = neuract(ind_AB,:);
                                    Echosen = neuract_OE(:,5) == -1; 
                                    averCritFR_OE = averCritFR(ind_AB,:);
                                    FR_in_OE = allneuroFR_icell(ind_AB,:);
                                end 
                                hits_chosenE = neuract(ind_chosenE,1);
                                hits_chosenO = neuract(ind_chosenO,1);
                                neuract_chosenE = actmake_TT(celldata, psyphydata, hits_chosenE, timewindow_new);
                                neuract_chosenO = actmake_TT(celldata, psyphydata, hits_chosenO, timewindow_new);                              
                                Echosen_ROC = [ones(size(neuract_chosenE)); zeros(size(neuract_chosenO))];
                                fr_act = [neuract_chosenE; neuract_chosenO];
                                %
                                trialset = {};
                                trialset(Echosen_ROC==1) = {'Eeasy'};
                                trialset(Echosen_ROC==0) = {'Oeasy'};
                                [~,~,~,auc] = perfcurve(trialset,fr_act,'Eeasy','nboot',100,'alpha',.05);
                                %
                                eval(['CP_in_OE_criTW_easy.',critTWname,'(:,ncells_ana) = auc(1);'])
                                eval(['CPpval_in_OE_criTW_easy.',critTWname,'(:,ncells_ana) = prod(auc(2:3)-.5)>0;'])
                                %                           
                                FR_in_OE_movTW_easy.Echosen(:,ncells_ana) = nanmean(FR_in_OE(Echosen==1,:));
                                FR_in_OE_movTW_easy.Ochosen(:,ncells_ana) = nanmean(FR_in_OE(Echosen==0,:));
                                %
                                % EO trials easy
                                if SOslope_icell == 1 % CJA
                                    ind_AB = orderABBA == 1 & ind_easy;
                                    ind_chosenE = ind_AB & neuract(:,5) ==  1;  
                                    ind_chosenO = ind_AB & neuract(:,5) == -1;  
                                    %
                                    neuract_EO = neuract(ind_AB,:);
                                    Echosen = neuract_EO(:,5) == 1;    
                                    averCritFR_EO = averCritFR(ind_AB,:);
                                    FR_in_EO = allneuroFR_icell(ind_AB,:);
                                elseif SOslope_icell == -1 %CJB
                                    ind_BA = orderABBA ==-1  & ind_easy;
                                    ind_chosenE = ind_BA & neuract(:,5) == -1;  
                                    ind_chosenO = ind_BA & neuract(:,5) ==  1; 
                                    %
                                    neuract_EO = neuract(ind_BA,:);
                                    Echosen = neuract_EO(:,5) == -1; 
                                    averCritFR_EO = averCritFR(ind_BA,:);
                                    FR_in_EO = allneuroFR_icell(ind_BA,:);
                                end 
                                hits_chosenE = neuract(ind_chosenE,1);
                                hits_chosenO = neuract(ind_chosenO,1);
                                neuract_chosenE = actmake_TT(celldata, psyphydata, hits_chosenE, timewindow_new);
                                neuract_chosenO = actmake_TT(celldata, psyphydata, hits_chosenO, timewindow_new);                              
                                Echosen_ROC = [ones(size(neuract_chosenE)); zeros(size(neuract_chosenO))];
                                fr_act = [neuract_chosenE; neuract_chosenO];
                                %
                                trialset = {};
                                trialset(Echosen_ROC==1) = {'Eeasy'};
                                trialset(Echosen_ROC==0) = {'Oeasy'};
                                [~,~,~,auc] = perfcurve(trialset,fr_act,'Eeasy','nboot',100,'alpha',.05);
                                %
                                eval(['CP_in_EO_criTW_easy.',critTWname,'(:,ncells_ana) = auc(1);'])
                                eval(['CPpval_in_EO_criTW_easy.',critTWname,'(:,ncells_ana) = prod(auc(2:3)-.5)>0;'])
                                %
                                FR_in_EO_movTW_easy.Echosen(:,ncells_ana) = nanmean(FR_in_EO(Echosen==1,:));
                                FR_in_EO_movTW_easy.Ochosen(:,ncells_ana) = nanmean(FR_in_EO(Echosen==0,:));
                                
                                
                            catch                               
                                eval(['CP_in_OE_criTW_hard.',critTWname,'(:,ncells_ana) = nan;'])
                                eval(['CPpval_in_OE_criTW_hard.',critTWname,'(:,ncells_ana) = 0;'])
                                %
                                eval(['CP_in_EO_criTW_hard.',critTWname,'(:,ncells_ana) = nan;'])
                                eval(['CPpval_in_EO_criTW_hard.',critTWname,'(:,ncells_ana) = 0;'])
                                %
                                eval(['CP_in_OE_criTW_easy.',critTWname,'(:,ncells_ana) = nan;'])
                                eval(['CPpval_in_OE_criTW_easy.',critTWname,'(:,ncells_ana) = 0;'])
                                %
                                eval(['CP_in_EO_criTW_easy.',critTWname,'(:,ncells_ana) = nan;'])
                                eval(['CPpval_in_EO_criTW_easy.',critTWname,'(:,ncells_ana) = 0;'])
                                %
                                FR_in_OE_movTW_hard.Echosen(:,ncells_ana) = nan([1, length(timepoints)]);
                                FR_in_OE_movTW_hard.Ochosen(:,ncells_ana) = nan([1, length(timepoints)]);
                                %
                                FR_in_EO_movTW_hard.Echosen(:,ncells_ana) = nan([1, length(timepoints)]);
                                FR_in_EO_movTW_hard.Ochosen(:,ncells_ana) = nan([1, length(timepoints)]);
                                %
                                FR_in_OE_movTW_easy.Echosen(:,ncells_ana) = nan([1, length(timepoints)]);
                                FR_in_OE_movTW_easy.Ochosen(:,ncells_ana) = nan([1, length(timepoints)]);
                                %
                                FR_in_EO_movTW_easy.Echosen(:,ncells_ana) = nan([1, length(timepoints)]);
                                FR_in_EO_movTW_easy.Ochosen(:,ncells_ana) = nan([1, length(timepoints)]);
                                %
                                eval(['CP_in_OEEO_criTW_hard.',critTWname,'(:,ncells_ana) = nan;'])
                                eval(['CPpval_in_OEEO_criTW_hard.',critTWname,'(:,ncells_ana) = 0;'])
                                %
                                eval(['CP_in_OEEO_criTW_easy.',critTWname,'(:,ncells_ana) = nan;'])
                                eval(['CPpval_in_OEEO_criTW_easy.',critTWname,'(:,ncells_ana) = 0;'])
                            end
                            
                        end % for icriTW
                    end % isequal JCorSO                   
                end %for icell 
                
                % save results
                filename = ['OEEOchoiceprop_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,dobhvlogit];   
                eval(['save ',filename ' monkeyname_list cellnames_list tuningslopes_list ncells_ana rho_all orderbias_all steepness_all timepoints corr_VOandVE '...
                                       ' CP_in_OE_movTW CP_in_EO_movTW CPpval_in_OE_movTW CPpval_in_EO_movTW '...
                                       ' CP_in_OE_criTW CP_in_EO_criTW CPpval_in_OE_criTW CPpval_in_EO_criTW '...
                                       ' CP_in_OE_criTW_noforced CP_in_EO_criTW_noforced CPpval_in_OE_criTW_noforced CPpval_in_EO_criTW_noforced '...
                                       ' CP_in_OE_criTW_hard CP_in_EO_criTW_hard CPpval_in_OE_criTW_hard CPpval_in_EO_criTW_hard '...
                                       ' CP_in_OE_criTW_easy CP_in_EO_criTW_easy CPpval_in_OE_criTW_easy CPpval_in_EO_criTW_easy '...
                                       ' FR_in_OE_movTW FR_in_EO_movTW FR_in_OE_movTW_noforced FR_in_EO_movTW_noforced '...
                                       ' FR_in_OE_movTW_hard FR_in_EO_movTW_hard FR_in_OE_movTW_easy FR_in_EO_movTW_easy '...
                                       ' CP_in_OEEO_criTW CPpval_in_OEEO_criTW CP_in_OEEO_criTW_noforced CPpval_in_OEEO_criTW_noforced '...
                                       ' CP_in_OEEO_criTW_hard CPpval_in_OEEO_criTW_hard CP_in_OEEO_criTW_easy CPpval_in_OEEO_criTW_easy '... 
                                       ])                    
               
            end % try catch
            
            
            % % % 
            % remove outlier - based on steepness
            ind_G = ismember(monkeyname_list,'G');
            ind_J = ismember(monkeyname_list,'J');
            ind_goodJC = (steepness_all.JC > steepoutlier_JC_G(1)  & steepness_all.JC < steepoutlier_JC_G(2)) & ind_G;
            ind_goodSO = (steepness_all.SO > steepoutlier_SO_G(1)  & steepness_all.SO < steepoutlier_SO_G(2)) & ind_G;
            % ind_noneout_G = ind_goodJC & ind_goodSO;
            ind_noneout_G = ind_goodSO;
            ind_goodJC = (steepness_all.JC > steepoutlier_JC_J(1)  & steepness_all.JC < steepoutlier_JC_J(2)) & ind_J;
            ind_goodSO = (steepness_all.SO > steepoutlier_SO_J(1)  & steepness_all.SO < steepoutlier_SO_J(2)) & ind_J;
            % ind_noneout_J = ind_goodJC & ind_goodSO; 
            ind_noneout_J = ind_goodSO; 
            ind_noneout = ind_noneout_G | ind_noneout_J;
            if ~dosteepout
                ind_noneout = logical(ones(size(ind_noneout)));
            end
            % % %
                     
            %
            %
            % PLOT - moving time window: average firing rate
            if 1
                if isequal(JCorSO,'SO')
                    %
                    CPtype = ''; % '', '_noforced', '_hard', '_easy'
                    CPtypename = 'hard trials'; % 'all trials', 'no forced choice', 'hard trials', 'easy trials'                   
                    %
                    if isequal(alignment,'off')
                        ind_time = timepoints>=-450 & timepoints<=2550;
                        eventpoints = [0, 500, 1000, 1500, 2000];
                        eventnames = {'offer1 on','offer1 off','offer2 on','offer2 off','target on'};
                    elseif isequal(alignment,'cho')
                        ind_time = timepoints>=-750 & timepoints<=750;
                        eventpoints = [0];
                        eventnames = {'outcome'};
                    end
                    XXX = timepoints(ind_time)';
                    eval(['YYY1 = nanmean(FR_in_EO_movTW',CPtype,'.Echosen(ind_time,ind_noneout)'');'])
                    eval(['YYY2 = nanmean(FR_in_EO_movTW',CPtype,'.Ochosen(ind_time,ind_noneout)'');'])
                    eval(['YYY3 = nanmean(FR_in_OE_movTW',CPtype,'.Echosen(ind_time,ind_noneout)'');'])
                    eval(['YYY4 = nanmean(FR_in_OE_movTW',CPtype,'.Ochosen(ind_time,ind_noneout)'');'])
                    XX = [min(XXX),max(XXX)];
                    % YY = [floor(min([YYY1,YYY2,YYY3,YYY4])), ceil(max([YYY1,YYY2,YYY3,YYY4]))];
                    YY = [4,9];
                    %
                    figure;
                    if isequal(alignment,'off')
                        set(gcf,'position',[110 65 1150 950], 'PaperPositionMode','auto')
                    elseif isequal(alignment,'cho')
                        set(gcf,'position',[110 65 650 950], 'PaperPositionMode','auto')
                    end
                    subplot(3,1,1);
                    hold on;
                    %
                    plot(XXX,YYY1,'Color',[0.0 0.0 0.6],'LineWidth',3);
                    plot(XXX,YYY2,'Color',[0.6 0.0 0.0],'LineWidth',3);
                    plot(XXX,YYY3,'Color',[0.0 1.0 1.0],'LineWidth',3);
                    plot(XXX,YYY4,'Color',[1.0 0.4 0.6],'LineWidth',3);
                    %
                    nevents = length(eventpoints);
                    for ievent = 1:nevents
                        eventpoint = [eventpoints(ievent),eventpoints(ievent)];
                        plot(eventpoint, YY, 'k--','LineWidth', 0.5); 
                        % text(eventpoints(ievent),YY(1)+(sum(YY)/25),eventnames{ievent},'fontsize',14);
                    end
                    set(gca,'XTick',eventpoints,'XTickLabel',eventnames);
                    %
                    if isequal(alignment,'off')
%                         twinnames_disp = {'before offers', 'late offer2','early wait','late wait','after wait'};
%                         twinranges = [-250 0; 1250 1500; 1500 1750; 1750 2000; 2000 2250];
%                         clorcodes = [0 0 0; 0 0 0; 0.35 0.35 0.35; 0.6 0.6 0.6; 0.8 0.8 0.8];
                        twinnames_disp = {'before offers', 'late offer2','early wait'};
                        twinranges = [-250 0; 1250 1500; 1500 1750];
                        clorcodes = [0 0 0; 0 0 0; 0.35 0.35 0.35];
                        for itwin = 1:length(twinnames_disp)
                            twinname_disp = twinnames_disp{itwin};
                            twinrange = twinranges(itwin,:);
                            plot(twinrange,[YY(1), YY(1)],'-','LineWidth',4,'color',clorcodes(itwin,:));
                            text(twinrange(1)+10,YY(1)+.25,{twinname_disp},'fontsize',12);
                        end
                    elseif isequal(alignment,'cho')
%                         twinnames_disp = {'early before outcome','before outcome','after outcome','later after outcome'};
%                         twinranges = [-500 -250; -250 0; 0 250; 250 500];
%                         clorcodes = [0 0 0; 0.35 0.35 0.35; 0.6 0.6 0.6; 0.8 0.8 0.8];
                        twinnames_disp = {'before outcome'};
                        twinranges = [-250 0];
                        clorcodes = [0.35 0.35 0.35];
                        for itwin = 1:length(twinnames_disp)
                            twinname_disp = twinnames_disp{itwin};
                            twinrange = twinranges(itwin,:);
                            plot(twinrange,[YY(1),YY(1)],'-','LineWidth',4,'color',clorcodes(itwin,:));
                            text(twinrange(1)+10,YY(1)+.25,{twinname_disp},'fontsize',12);
                        end
                    end
                    %
                    legend({'EO trial; E chosen'; 'EO trial; O chosen'; 'OE trial; E chosen'; 'OE trial; O chosen'})
                    xlabel('time(ms)');
                    ylabel('firing rate (sp/s)');
                    axis([XX YY]);
                    box off
                    set(gca,'FontSize',17);
   
                    % 
                    % moving time window ROC
                    subplot(3,1,2);
                    hold on;
                    XXX = timepoints(ind_time)';
                    YYY = nanmean(CP_in_OE_movTW(ind_time,ind_noneout)');
                    p_ttest = [];
                    for itp = 1:size(CP_in_OE_movTW,1)
                        [~, p_ttest(1,itp)] = ttest(CP_in_OE_movTW(itp,:),0.5);
                    end
                    XXXstars = XXX(p_ttest(ind_time)<0.005);
                    XX = [min(XXX),max(XXX)];
                    YY = [0.44 0.56];
                    %
                    plot(XXX,YYY,'Color',[0.0 0.0 0.0],'LineWidth',3);
                    plot(XXX,0.5*ones(size(XXX)),'k--','LineWidth',1.5);
                    plot(XXXstars,0.55*ones(size(XXXstars)),'r.');
                    %
                    nevents = length(eventpoints);
                    for ievent = 1:nevents
                        eventpoint = [eventpoints(ievent),eventpoints(ievent)];
                        plot(eventpoint, YY, 'k--','LineWidth', 0.5);
                        % text(eventpoints(ievent),YY(1)+(sum(YY)/25),eventnames{ievent},'fontsize',14);
                    end
                    set(gca,'XTick',eventpoints,'XTickLabel',eventnames);
                    %
                    legend({'ROC of OE trials'},'Location','southeast')
                    xlabel('time(ms)');
                    ylabel('choice probability');
                    axis([XX YY]);
                    box off
                    set(gca,'FontSize',17);      
                    
                     % 
                    % moving time window ROC (CP) V.S. preference bias index
                    subplot(3,1,3);
                    hold on;
                    XXX = timepoints(ind_time)';
                    YYY  = [];
                    p_YYY = [];
                    norm_delta_rho = (rho_all.SO - rho_all.JC)./(rho_all.SO + rho_all.JC);
                    norm_delta_rho = norm_delta_rho(ind_noneout);
                    CP_tgt = CP_in_OE_movTW(ind_time,ind_noneout);
                    for itw = 1:sum(ind_time)
                        CP_tgt_itw = CP_tgt(itw,:)';
                        [YYY(1,itw), p_YYY(1,itw)] = corr(CP_tgt_itw(~isnan(CP_tgt_itw)),norm_delta_rho(~isnan(CP_tgt_itw))','type','Spearman');
                    end
                    %
                    XXXstars = XXX(p_YYY<0.05);
                    XX = [min(XXX),max(XXX)];
                    YY = [-0.3 0.3];
                    %
                    plot(XXX,YYY,'Color',[0.0 0.0 0.0],'LineWidth',3);
                    plot(XXX,0.00*ones(size(XXX)),'k--','LineWidth',1.5);
                    plot(XXXstars,0.28*ones(size(XXXstars)),'r.');
                    %
                    nevents = length(eventpoints);
                    for ievent = 1:nevents
                        eventpoint = [eventpoints(ievent),eventpoints(ievent)];
                        plot(eventpoint, YY, 'k--','LineWidth', 0.5);
                        % text(eventpoints(ievent),YY(1)+(sum(YY)/25),eventnames{ievent},'fontsize',14);
                    end
                    set(gca,'XTick',eventpoints,'XTickLabel',eventnames);
                    %
                    legend({'corr between CP and PBI'},'Location','southeast')
                    xlabel('time(ms)');
                    ylabel('corr');
                    axis([XX YY]);
                    box off
                    set(gca,'FontSize',17);      
                    
                    %
                    %
                    axes('position',[.02 .95 .2 .05]);
                    text(0,0,{['chosen juice cells'];['Task 2 in Two Task'];CPtypename},'fontsize',10);
                    axis off
                    
                end % isequal JCorSO
            end
            
            
            %
            %
            % PLOT - critical time window
            if 1
                if isequal(JCorSO,'SO')
                    %
                    CPtype = '_hard'; % '', '_noforced', '_hard', '_easy'
                    CPtypename = 'hard trials'; % 'all trials', 'no forced choice', 'hard trials', 'easy trials'
                    % 
                    if isequal(alignment,'off')
                        twinnames = {'beforeofferon','lateoffer2','earlywait','latewait','afterwait'};
                        twinnames_disp = {'before offers','late offer2','early wait','late wait','after wait'};
                        TWperiods = {'-250~0ms before offers','250-500ms after offer2 onset', '0-250ms in wait', '250-500ms in wait', '0-250ms after wait'};
                    elseif isequal(alignment,'cho')
                        twinnames = {'earlyprejuice', 'prejuice', 'postjuice','latepostjuice'};
                        twinnames_disp = {'early before outcome','before outcome','after outcome','later after outcome'};
                        TWperiods = {'-500~-250ms before outcome', '-250~0ms before outcome', '0-250ms after outcome', '250-500ms after outcome'};
                    end
                    %
                    for itwin = 1:length(twinnames)
                        twinname = twinnames{itwin};
                        eval(['CP_in_merged_criTW',CPtype,'.',twinname,' = [CP_in_OE_criTW',CPtype,'.',twinname,', CP_in_EO_criTW',CPtype,'.',twinname,'];'])
                        eval(['CPpval_in_merged_criTW',CPtype,'.',twinname,' = [CPpval_in_OE_criTW',CPtype,'.',twinname,', CPpval_in_EO_criTW',CPtype,'.',twinname,'];'])                   
                    end
                               
                    for itwin = 1:length(twinnames)
                        twinname = twinnames{itwin};
                        twinname_disp = twinnames_disp{itwin};
                        TWperiod = TWperiods{itwin};
                        %
                        figure;
                        set(gcf,'position',[110 65 1350 1550], 'PaperPositionMode','auto')
                        %
%                         histtypes = {'OE', 'EO', 'merged'};
                        histtypes = {'OE'};
                        %
                        % xxxtypes = {'rho_all.SO','steepness_all.SO','orderbias_all.SO'};
                        % xlabelnames = {'\rho Task 2', '\eta Task 2', '\epsilon Task 2'};
                        norm_delta_rho = (rho_all.SO - rho_all.JC)./(rho_all.SO + rho_all.JC);
                        % norm_delta_rho = (rho_all.SO - rho_all.JC);
                        delta_steepness = steepness_all.SO - steepness_all.JC;
%                         xxxtypes = {'norm_delta_rho','delta_steepness','orderbias_all.SO'};
%                         xlabelnames = {'normalized \Delta\rho', '\Delta\eta (Task 2-Task 1)', '\epsilon Task 2'};
%                         yyytypes = {'OE','EO','merged'};
                        xxxtypes = {'norm_delta_rho'};
                        xlabelnames = {'normalized \Delta\rho'};
                        yyytypes = {'OE'};
                        
                        for ihist = 1:length(histtypes)
                            histtype = histtypes{ihist};
                            % 
                            % eval(['ind_good = logical(ones(size(CP_in_',histtype,'_criTW',CPtype,'.',twinname,')));']);
                            % eval(['ind_good = CPpval_in_',histtype,'_criTW',CPtype,'.',twinname,';']);  
                            % eval(['ind_good = CPpval_in_',histtype,'_criTW',CPtype,'.',twinname,' & CP_in_',histtype,'_criTW',CPtype,'.',twinname '>0.5;']);  
                            ind_good = ind_noneout;
                            %
                            subplot(length(xxxtypes)+1,length(histtypes),ihist);
                            eval(['XXX = [CP_in_',histtype,'_criTW',CPtype,'.',twinname,'(ind_good)];'])
                            %
                            ind_nan = isnan(XXX);
                            XXX = XXX(~ind_nan);
                            %
                            XX = [0.05, 0.95]; 
                            YY = [0, 50];
                            edges = [XX(1):0.05:XX(2)];
                            hold on; 
                            plot([0.5 0.5],YY, '--','Color', [0.5 0.5 0.5],'LineWidth',1); hold on;
                            histogram(XXX,edges,'FaceColor',[.4 .4 .4]); hold on;                   
                            [~, pp_ttest] = ttest(XXX,0.5);
                            [ pp_wil,  ~] = signrank(XXX,0.5);
                            title([histtype,' trials in ',twinname_disp])
                            xlabel('choice probability');
                            ylabel('cell #');
                            text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
                                {['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                                 ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                                 ['mean = ',num2str(nanmean(XXX),'%.4f')];... 
                                 ['std = ',num2str(nanstd(XXX),'%.4f')];... 
                                 },'fontsize',12);
                            box off
                            axis([XX YY])
                            axis square    
                            set(gca,'FontSize',16); 
                        end
                    
                        %
                        for ixtype = 1:length(xxxtypes)
                            
                            for iytype = 1:length(yyytypes)
                                yyytype = yyytypes{iytype};
                                %
                                % eval(['ind_good = logical(ones(size(CP_in_',yyytype,'_criTW',CPtype,'.',twinname,')));']);
                                % eval(['ind_good = CPpval_in_',yyytype,'_criTW',CPtype,'.',twinname,';']);
                                % eval(['ind_good = CPpval_in_',yyytype,'_criTW',CPtype,'.',twinname,' & CP_in_',yyytype,'_criTW',CPtype,'.',twinname '>0.5;']);  
                                ind_good = ind_noneout;
                                %
                                subplot(length(xxxtypes)+1,length(histtypes),iytype+length(histtypes)*ixtype)
                                %
                                xxxtype = xxxtypes{ixtype};
                                xlabelname = xlabelnames{ixtype};
                                eval(['XXX = ',xxxtype,';'])                           
                                if isequal(yyytype,'merged')
                                    XXX = [XXX,XXX];
                                end
                                XXX = XXX(ind_good);                           
                                eval(['YYY = [CP_in_',yyytype,'_criTW',CPtype,'.',twinname,'(ind_good)];'])
                                %
                                ind_nan = isnan(XXX) | isnan(YYY);
                                XXX = XXX(~ind_nan);
                                YYY = YYY(~ind_nan);
                                %
                                % remove outlier
%                                 kout = 1.5;    
%                                 quantiles_XXX = quantile(XXX,[0.25 0.5 0.75]);
%                                 IQR_XXX = quantiles_XXX(3) - quantiles_XXX(1);
%                                 outlier_XXX = [quantiles_XXX(1)-kout*IQR_XXX, quantiles_XXX(3)+kout*IQR_XXX];
%                                 ind_OUT_XXX = XXX > outlier_XXX(2) | XXX < outlier_XXX(1);
%                                 XXX = XXX(~ind_OUT_XXX); 
%                                 YYY = YYY(~ind_OUT_XXX); 
                                
                                % K = 3;                            
                                % ind_OUT_XXX = (XXX>(mean(XXX)+std(XXX)*K)) | (XXX<(mean(XXX)-std(XXX)*K));  
                                % ind_OUT_XXX = ind_OUT_XXX | XXX < -15; % for steepess
                                % XXX = XXX(~ind_OUT_XXX); 
                                % YYY = YYY(~ind_OUT_XXX); 
                                %
                                aa = polyfit(XXX,YYY,1);
                                % XX = [floor(min(XXX)),ceil(max(XXX))];
                                XX = [min(XXX)-(max(XXX)-min(XXX))/10, max(XXX)+(max(XXX)-min(XXX))/10];
                                YY = [floor(min(YYY)),ceil(max(YYY))]; 
                                %
                                aa = polyfit(XXX,YYY,1);
                                Yfit = aa(1)*XX+aa(2);
                                %
                                hold on; plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
                                hold on; plot(XXX, YYY, 'ko','MarkerSize',10);
                                Sigma_ell = cov(XXX, YYY);
                                mu_ell(1) = nanmean(XXX);
                                mu_ell(2) = nanmean(YYY);       
                                hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                                title([yyytype,' trials in ',twinname_disp])
                                xlabel(xlabelname);
                                ylabel('choice probability');
                                [RR_Pea, pp_Pea] = corr(XXX', YYY', 'type', 'Pearson');
                                [RR_Spe, pp_Spe] = corr(XXX', YYY', 'type', 'Spearman');
                                if pp_Pea<0.05 | pp_Spe<0.05
                                    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
                                    {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
                                     ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
                                     ['N = ',num2str(length(XXX)),' responses']}, 'fontsize', 12, 'color', 'k');
                                else
                                    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
                                    {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
                                     ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
                                     ['N = ',num2str(length(XXX)),' responses']}, 'fontsize', 12, 'color', 'k');
                                end
                                box off
                                axis([XX YY])
                                axis square    
                                set(gca,'FontSize',16); 
                            
                            end % for iytype
                        end % for ixtype
                        
                        % 
                        axes('position',[.02 .95 .2 .05]);
                        text(0,0,{['chosen juice cells'];TWperiod;['Task 2 in Two Task'];CPtypename},'fontsize',10);
                        axis off
                    
                    end % for itwin                  
                end % isequal JCSO
                
                %
                % behavioral measurement comparison
                %
                figure;
                set(gcf,'position',[110 65 1350 1550], 'PaperPositionMode','auto')
                %  
                norm_delta_rho = (rho_all.SO - rho_all.JC)./(rho_all.SO + rho_all.JC);
                delta_steepness = steepness_all.SO - steepness_all.JC;                
                xxxtypes = {'orderbias_all.SO', 'orderbias_all.SO', 'norm_delta_rho'};  
                yyytypes = {'norm_delta_rho',   'delta_steepness',  'delta_steepness'};
                xlabelnames = {'\epsilon Task 2',       '\epsilon Task 2',            'normalized \Delta\rho'};
                ylabelnames = {'normalized \Delta\rho', '\Delta\eta (Task 2-Task 1)', '\Delta\eta (Task 2-Task 1)'};
                ntypes = 3;
                %
                for itype = 1:ntypes                  
                    %
                    subplot(1,ntypes,itype)
                    %
                    xxxtype = xxxtypes{itype};
                    xlabelname = xlabelnames{itype};
                    eval(['XXX = ',xxxtype,';'])    
                    %
                    yyytype = yyytypes{itype};
                    ylabelname = ylabelnames{itype};
                    eval(['YYY = ',yyytype,';'])   
                    %
                    % remove outlier
                    kout = 1.5;    
                    quantiles_XXX = quantile(XXX,[0.25 0.5 0.75]);
                    IQR_XXX = quantiles_XXX(3) - quantiles_XXX(1);
                    outlier_XXX = [quantiles_XXX(1)-kout*IQR_XXX, quantiles_XXX(3)+kout*IQR_XXX];
                    ind_OUT_XXX = XXX > outlier_XXX(2) | XXX < outlier_XXX(1);
                    quantiles_YYY = quantile(YYY,[0.25 0.5 0.75]);
                    IQR_YYY = quantiles_YYY(3) - quantiles_YYY(1);
                    outlier_YYY = [quantiles_YYY(1)-kout*IQR_YYY, quantiles_YYY(3)+kout*IQR_YYY];
                    ind_OUT_YYY = YYY > outlier_YYY(2) | YYY < outlier_YYY(1);
                    ind_OUT = ind_OUT_XXX | ind_OUT_YYY;
                    XXX = XXX(~ind_OUT); 
                    YYY = YYY(~ind_OUT); 
                    %
                    % K = 3;                            
                    % ind_OUT_XXX = (XXX>(mean(XXX)+std(XXX)*K)) | (XXX<(mean(XXX)-std(XXX)*K));  
                    % ind_OUT_XXX = ind_OUT_XXX | XXX < -15; % for steepess
                    % ind_OUT_YYY = (YYY>(mean(YYY)+std(YYY)*K)) | (YYY<(mean(YYY)-std(YYY)*K));  
                    % ind_OUT_YYY = ind_OUT_YYY | YYY < -15; % for steepess
                    % ind_OUT = ind_OUT_XXX | ind_OUT_YYY;
                    % XXX = XXX(~ind_OUT); 
                    % YYY = YYY(~ind_OUT); 
                    %
                    XX = [min(XXX)-(max(XXX)-min(XXX))/10, max(XXX)+(max(XXX)-min(XXX))/10];
                    YY = [min(YYY)-(max(YYY)-min(YYY))/10, max(YYY)+(max(YYY)-min(YYY))/10];
                    %
                    aa = deming(XXX',YYY');
                    Yfit = aa(2)*XX+aa(1);
                    %
                    hold on; plot(XX,Yfit,'-','LineWidth',2, 'Color', [0.4 0.4 0.4]);
                    hold on; plot(XXX, YYY, 'ko','MarkerSize',6);
                    Sigma_ell = cov(XXX, YYY);
                    mu_ell(1) = nanmean(XXX);
                    mu_ell(2) = nanmean(YYY);       
                    hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);                    
                    xlabel(xlabelname);
                    ylabel(ylabelname);
                    [RR_Pea, pp_Pea] = corr(XXX', YYY', 'type', 'Pearson');
                    [RR_Spe, pp_Spe] = corr(XXX', YYY', 'type', 'Spearman');
                    %
                    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
                    {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
                     ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
                     ['N = ',num2str(length(XXX)),' responses']}, 'fontsize', 11, 'color', 'k');
                    %   
                    box off
                    axis([XX YY])
                    axis square    
                    set(gca,'FontSize',14); 
                end   
            end
            
   
                        
            
        end %for ialign
    end %for iJCSO  
end %for iclass



%%
function plotErrorEllipse(mu_ell, Sigma_ell, p_ell)
    s = -2 * log(1 - p_ell);
    [V, D] = eig(Sigma_ell * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    plot(a(1, :) + mu_ell(1), a(2, :) + mu_ell(2), 'LineWidth',1.5, 'Color',[0.4 0.4 0.4] );
end
