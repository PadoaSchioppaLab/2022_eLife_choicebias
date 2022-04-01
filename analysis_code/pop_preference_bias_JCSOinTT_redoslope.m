% pop_preference_bias_JCSOinTT.m
%
% This script analyzes the neural correlates for the preference bias in JC and SO.
% This script focuses on postoffer time window for JC and postoffer1 and postoffer2 time window for SO
% offer value cells  - EO (offer1 TW) v.s. OE (offer2 TW): slope and FR
% 
% This script uses plot_pop_preference_bias_JCSOinTT.m to plot 
%

%
% author:     Feb 2021; WS


close all
clearvars

brainarea = 'OFC'; % 'DLPFC', 'VLPFC' 'OFC'
monkey_ana = 'both'; % 'Gervinho', 'Juan', 'both'

% % % 
JCSOclassifySep = 0; % if 1, load neurons that is classified seperately by JC or SO trials: 0-JCSOcellist

slopepooled = 1;     % if 1, pool cell postive and negative

atleast_nntrials = 2;

dobhvlogit = '_logit'; % '': probit; '_logit'; '_neworderbias'

saveplots = 0;
if saveplots
if slopepooled
    figuresave = ['C:\Experiments\TwoTasks\Analysis\Analysis_OFC\pop_preference_orderbias_JCSOinTT_redoslope_slopepooled_',monkey_ana];
elseif ~slopepooled
    figuresave = ['C:\Experiments\TwoTasks\Analysis\Analysis_OFC\pop_preference_orderbias_JCSOinTT_redoslope_',monkey_ana];    
end
eval(['!del ',figuresave,'.ps; !del ',figuresave,'.pdf']);
end

% % % 
% generate contingency table with TTcellist_OFC_xx.m
differentRho = '_differRho';    % '_differRho' or ''; do analysis based on the same rho in JC and SO or not 
twoTWinJCandSO = '';  % '' or ''; do fewer TW: 2 for JC and 2 for SO
doJCseq = '';          % '_JCseq' or ''; do sequential JC parameters or not
nntrials = '_2nntrials'; % '_2nntrials' or '_3nntrials'
anovasetup = '_01pvallessTW_95nonzero'; % '_05pval_95nonzero'; '_01pvallessTW_95nonzero';
filename = ['TTcellist',doJCseq,'_',brainarea,'_both',differentRho,nntrials,anovasetup];
load(filename);

% % %
% example neurons
examplecells = {'J190806a41'}; %,... % OV: 'G190112a43', 'J191111b33', 'G181209b31', 'J190806a41', 'G181025b11', 'J190721c32', 'G181206b42'            
                %'G190103b41',... % CJ: 'G190103b41',
                % 'J191122a43'};   % CV: 'G181018a11', 'G181208a33', 'G181208a31'
doexamplecell = 0;

% % % % % % 
% remove outlier based on steepness across all sessions (keep consistent with behavioral measurement)
removesessions = 0; % remove based on dynamic range and saturation
dosteepout = 1; 
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


% % % 
cellnames = {JCcellist.infocell.cellname};
monkeynames = [];
for icellname = 1:length(cellnames)
    monkeynames{icellname,1} = cellnames{icellname}(1);
end
if isequal(monkey_ana,'Gervinho')
    ind_G = ismember(monkeynames,'G');
    JCcellist.infocell = JCcellist.infocell(ind_G);
    JCcellist.specs.ncells = sum(ind_G);
    SOcellist.infocell = SOcellist.infocell(ind_G);
    SOcellist.specs.ncells = sum(ind_G);
    JCSOcellist.infocell = JCSOcellist.infocell(ind_G);
    JCSOcellist.specs.ncells = sum(ind_G);
elseif isequal(monkey_ana,'Juan')
    ind_J = ismember(monkeynames,'J');
    JCcellist.infocell = JCcellist.infocell(ind_J);
    JCcellist.specs.ncells = sum(ind_J);
    SOcellist.infocell = SOcellist.infocell(ind_J);
    SOcellist.specs.ncells = sum(ind_J);
    JCSOcellist.infocell = JCSOcellist.infocell(ind_J);
    JCSOcellist.specs.ncells = sum(ind_J);
end

% initiation
allmonkeys    = [];
allcellnames  = [];
allFRs        = [];
allSteepness  = [];
allOrderbias  = [];
allsessrange  = [];
allrhos       = [];
%
allslopes     = [];
allnonzeros   = [];
allRsq        = [];
allintcepts   = [];
allmaxinters  = [];
allslopes_QAB = [];     % for CV cells
allintcepts_QAB = []; % for CV cells, postjuice TW 
allnonzeros_QAB = [];   % for CV cells


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
%
cellclassnames = {'OVA', 'OVB'};
cellclassnums = [1,2];
nclassnames = size(cellclassnames,2);
%
slopesigns = [1, -1];
slopesignnames = {'positive','negative'};
nslopesigns = size(slopesigns,2);
%

for iclassname = 1:nclassnames
    classname = cellclassnames{iclassname};
    %
    for islopesign = 1:nslopesigns
        slopesignname = slopesignnames{islopesign};
        %
        try
            exnovo
            filename = ['pop_preference_bias_JCSOinTT_redoslope',differentRho,'_',monkey_ana,dobhvlogit];            
            load(filename);   
            %
            eval(['cellnames_iclass = allcellnames.',classname,'.',slopesignname,';'])
            eval(['monkeyname_list = allmonkeys.',classname,'.',slopesignname,';'])
            %
            eval(['FR_inJC = allFRs.',classname,'.',slopesignname,'.JC;'])
            eval(['FR_inSO = allFRs.',classname,'.',slopesignname,'.SO;'])
            eval(['steepness_inJC = allSteepness.',classname,'.',slopesignname,'.JC;'])
            eval(['steepness_inSO = allSteepness.',classname,'.',slopesignname,'.SO;'])
            eval(['orderbias_inSO = allOrderbias.',classname,'.',slopesignname,'.SO;'])
            eval(['sessrange_inJC = allsessrange.',classname,'.',slopesignname,'.JC;'])
            eval(['sessrange_inSO = allsessrange.',classname,'.',slopesignname,'.SO;'])
            eval(['rho_inJC = allrhos.',classname,'.',slopesignname,'.JC;'])
            eval(['rho_inSO = allrhos.',classname,'.',slopesignname,'.SO;'])
            % 
            eval(['slope_inJC = allslopes.',classname,'.',slopesignname,'.JC;'])
            eval(['slope_inSO = allslopes.',classname,'.',slopesignname,'.SO;'])        
            eval(['Rsq_inJC = allRsq.',classname,'.',slopesignname,'.JC;'])
            eval(['Rsq_inSO = allRsq.',classname,'.',slopesignname,'.SO;'])
            eval(['nonzero_inJC = allnonzeros.',classname,'.',slopesignname,'.JC;'])
            eval(['nonzero_inSO = allnonzeros.',classname,'.',slopesignname,'.SO;'])            
            eval(['intcept_inJC = allintcepts.',classname,'.',slopesignname,'.JC;'])
            eval(['intcept_inSO = allintcepts.',classname,'.',slopesignname,'.SO;'])
            eval(['maxintc_inJC = allmaxinters.',classname,'.',slopesignname,'.JC;'])
            eval(['maxintc_inSO = allmaxinters.',classname,'.',slopesignname,'.SO;'])
            %
            eval(['slope_QAB_inJC = allslopes_QAB.',classname,'.',slopesignname,'.JC;'])
            eval(['intcept_QAB_inJC = allintcepts_QAB.',classname,'.',slopesignname,'.JC;'])
            eval(['nonzero_QAB_inJC = allnonzeros_QAB.',classname,'.',slopesignname,'.JC;']) 
            eval(['slope_QAB_inSO = allslopes_QAB.',classname,'.',slopesignname,'.SO;'])
            eval(['intcept_QAB_inSO = allintcepts_QAB.',classname,'.',slopesignname,'.SO;'])
            eval(['nonzero_QAB_inSO = allnonzeros_QAB.',classname,'.',slopesignname,'.SO;'])          
        catch
            
            ind_JCiclass = JCclasses==cellclassnums(iclassname) & JCslopes==slopesigns(islopesign);
            ind_SOiclass = SOclasses==cellclassnums(iclassname) & SOslopes==slopesigns(islopesign);
            ind_iclass = ind_JCiclass & ind_SOiclass;
            % ind_iclass = ind_SOiclass;
            %
            ncells = sum(ind_iclass);
            cellnames_iclass = cellnames(ind_iclass);
            cellnum_iclass = find(ind_iclass==1);

            % innitiation for each condition
            monkeyname_list = [];
            %
            steepness_inJC = []; steepness_inSO = [];
            rho_inJC = [];       rho_inSO = [];
            sessrange_inJC = []; sessrange_inSO = [];
            orderbias_inSO = [];
            %
            slope_inJC = [];     slope_inSO = [];
            Rsq_inJC = [];       Rsq_inSO = [];
            nonzero_inJC = [];   nonzero_inSO = [];
            FR_inJC = [];        FR_inSO = [];
            intcept_inJC = [];   intcept_inSO = [];
            maxintc_inJC = [];   maxintc_inSO = [];
            valrange_inJC = [];  valrange_inSO = [];
            %     
            slope_QAB_inSO = [];
            intcept_QAB_inSO = [];
            nonzero_QAB_inSO = [];
            slope_QAB_inJC = [];
            intcept_QAB_inJC = [];
            nonzero_QAB_inJC = [];
            
            for icell = 1:ncells           

                cellname = cellnames_iclass{icell};
                session = cellname(1:8);
                readsession_TT;
                % 
                monkeyname_list{icell} = cellname(1);
                    
                % % % 
                % load behavior results                              
                try
                    exnovo
                    filename = [dirroot,cellname,'_psyphycell'];
                    eval(['load ',filename])                                    
                    if isempty(differentRho)
                        rho_inJC(icell,:) = (psyphycell.sigmoidfit.JC{3}(1)+psyphycell.sigmoidfit.SO{3}(1))/2;
                        rho_inSO(icell,:) = (psyphycell.sigmoidfit.JC{3}(1)+psyphycell.sigmoidfit.SO{3}(1))/2;
                    elseif ~isempty(differentRho)
                        rho_inJC(icell,:) = psyphycell.sigmoidfit.JC{3}(1);                        
                        rho_inSO(icell,:) = psyphycell.sigmoidfit.SO{3}(1); 
                    end
                    steepness_inJC(icell,:) = psyphycell.sigmoidfit.JC{2}(2);
                    steepness_inSO(icell,:) = psyphycell.sigmoidfit.SO{2}(2);
                    orderbias_inSO(icell,:) = -2.*rho_inSO(icell,:).*psyphycell.sigmoidfit.SO{2}(3)./psyphycell.sigmoidfit.SO{2}(2);    
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
                        rho_inJC(icell,:) = (relvalue_SO + relvalue_JC)/2;
                        rho_inSO(icell,:) = (relvalue_SO + relvalue_JC)/2;
                    elseif ~isempty(differentRho)
                        rho_inJC(icell,:) = relvalue_JC;
                        rho_inSO(icell,:) = relvalue_SO;
                    end 
                    steepness_inJC(icell,:) = psyphycell.JC.NonChHyst.sigmoidfit.beta(2);
                    steepness_inSO(icell,:) = psyphycell.SO.NonChHyst.sigmoidfit.beta(2);  
                    if isequal(dobhvlogit, '_neworderbias')
                        rho_inAB = exp(-(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(2)));
                        rho_inBA = exp(-(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(2)));
                        orderbias_inSO(icell,:) = rho_inBA - rho_inAB;
                    else
                        orderbias_inSO(icell,:) = -2.*rho_inSO(icell,:).*psyphycell.SO.NonChHyst.sigmoidfit.beta(3)./psyphycell.SO.NonChHyst.sigmoidfit.beta(2); 
%                       orderbias_inSO(icell,:) = -2.*rho_inSO(icell,:).*psyphycell.SO.OrdChHyst.sigmoidfit.beta(3)./psyphycell.SO.OrdChHyst.sigmoidfit.beta(2); 
                    end
                    
                end                
                %
                
                % % % 
                filename = [dirroot,cellname,'_tuning'];
                eval(['load ',filename])
                
                % redo the linear regression
                % one time windows: JC: postoffer               
                twins = {'postoffer'};
                ntwins = length(twins);
                for iTW = 1:ntwins
                    twin = twins{iTW};
                    eval(['neuract_JC = tuning.JC.AB.neuract.bytrial.',twin,';'])
                    act_JC = neuract_JC(:,7);
                    OVA_JC = neuract_JC(:,2).*rho_inJC(icell);
                    OVB_JC = neuract_JC(:,3);
                    valrange_OVA_JC = nanmax(neuract_JC(:,2).*rho_inJC(icell));
                    valrange_OVB_JC = nanmax(neuract_JC(:,3));
                    CJ_JC = neuract_JC(:,4).*neuract_JC(:,5);  % -1: B
                    CV_JC = neuract_JC(:,2).*rho_inJC(icell);
                    CV_JC(CJ_JC==-1,:) = neuract_JC(CJ_JC==-1,3);                
                    valrange_CV_JC = nanmax(CV_JC) - nanmin(CV_JC);
                    sessrange_inJC(icell,iTW) = sqrt(valrange_OVA_JC*valrange_OVB_JC);
                    if     isequal(classname,'OVA'), valrange_inJC(icell,iTW) = valrange_OVA_JC;
                    elseif isequal(classname,'OVB'), valrange_inJC(icell,iTW) = valrange_OVB_JC;
                    elseif isequal(classname,'CV'), valrange_inJC(icell,iTW) = valrange_CV_JC;    
                    elseif isequal(classname,'CJ'), valrange_inJC(icell,iTW) = 1; end  
                    %
                    % average over trial type: [OVA, OVB, CJ]
                    % add other parameters into trial type table
                    % [OVA, OVB, CJ, #trials, propB, FR, std, relvalue]
                    newTrialTypes_all = [OVA_JC, OVB_JC, CJ_JC];
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
                        trialnum_chB = sum(CJ_JC(ind_newTriType)==-1);
                        TrialTypes_tbl(inewTriType,5) = trialnum_chB./trialnum;
                        %
                        FR_all = act_JC(ind_newTriType,:);
                        FR_itritype = nanmean(FR_all);
                        FR_std_itritype = nanstd(FR_all);
                        TrialTypes_tbl(inewTriType,6) = FR_itritype; 
                        TrialTypes_tbl(inewTriType,7) = FR_std_itritype;  
                        TrialTypes_tbl(inewTriType,8) = rho_inJC(icell);      
                    end
                    ind_badtrialtype = TrialTypes_tbl(:,4)<atleast_nntrials;
                    TrialTypes_tbl(ind_badtrialtype,:) = [];                                                     
                    %
                    % define independent variables (use same name as cell type)                       
                    OVA = TrialTypes_tbl(:,1);
                    OVB = TrialTypes_tbl(:,2);
                    CJ = TrialTypes_tbl(:,3);   % -1 B; 1 A    
                    CV = OVA;
                    CV(CJ ==-1) = OVB(CJ ==-1);  
                    % fitting
                    eval(['XXX = ',classname,';'])
                    YYY = TrialTypes_tbl(:,6);
                    g = fittype('m*x+q','coeff',{'m','q'});
                    [linfit,goodness] = fit(XXX,YYY,g);
                    %
                    FR_inJC(icell,iTW) = nanmean(act_JC);              
                    Rsq_inJC(icell,iTW) = goodness.rsquare;
                    slope_inJC(icell,iTW) = linfit.m;
                    intcept_inJC(icell,iTW) = linfit.q;
                    interval = confint(linfit,0.90);
                    if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                    nonzero_inJC(icell,iTW) = nonzero;                
                    maxintc_inJC(icell,iTW) = slope_inJC(icell,iTW).*valrange_inJC(icell,iTW)+intcept_inJC(icell,iTW);
                    %
                    % for offer value cells, fit quantity of A or B
                    if ismember(classname,{'OVA','OVB'})
                        OVA = TrialTypes_tbl(:,1)./rho_inJC(icell);
                        OVB = TrialTypes_tbl(:,2);
                        eval(['XXX = ',classname,';'])
                        YYY = TrialTypes_tbl(:,6);
                        g = fittype('m*x+q','coeff',{'m','q'});
                        [linfit,goodness] = fit(XXX,YYY,g);
                        %
                        slope_QAB_inJC(icell,iTW) = linfit.m;
                        intcept_QAB_inJC(icell,iTW) = linfit.q;
                        interval = confint(linfit,0.90);
                        if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                        nonzero_QAB_inJC(icell,iTW) = nonzero; 
                    else
                        slope_QAB_inJC = [];
                        nonzero_QAB_inJC = [];
                        intcept_QAB_inJC = [];
                    end
                end
                
                
                % %
                % redo the linear regression
                % two time windows: SO: postoffer1 postoffer2
                twins = {'postoffer1','postoffer2'};
                ntwins = length(twins);
                for iTW = 1:ntwins
                    twin = twins{iTW};
                    eval(['neuract_SO = tuning.SO.ABA.neuract.bytrial.',twin,';'])
                    act_SO = neuract_SO(:,7);
                    OVA_SO = neuract_SO(:,2).*rho_inSO(icell);
                    OVB_SO = neuract_SO(:,3);
                    valrange_OVA_SO = nanmax(neuract_SO(:,2).*rho_inSO(icell));
                    valrange_OVB_SO = nanmax(neuract_SO(:,3));
                    CJ_SO = neuract_SO(:,5);  % -1: B
                    ord_SO = neuract_SO(:,6); % -1: BA
                    CO_SO = CJ_SO.*ord_SO;    % -1: chosen 2nd
                    CV_SO = neuract_SO(:,2).*rho_inSO(icell);
                    CV_SO(CJ_SO==-1,:) = neuract_SO(CJ_SO==-1,3);
                    valrange_CV_SO = nanmax(CV_SO) - nanmin(CV_SO);
                    sessrange_inSO(icell,iTW) = sqrt(valrange_OVA_SO*valrange_OVB_SO);               
                    if     isequal(classname,'OVA'), valrange_inSO(icell,iTW) = valrange_OVA_SO;
                    elseif isequal(classname,'OVB'), valrange_inSO(icell,iTW) = valrange_OVB_SO;
                    elseif isequal(classname,'CV'), valrange_inSO(icell,iTW) = valrange_CV_SO; 
                    elseif isequal(classname,'CJ'), valrange_inSO(icell,iTW) = 1; end                
                    %
                    % average over trial type: [OVA, OVB, ord, CJ]
                    % add other parameters into trial type table
                    % [OVA, OVB, ord, CJ, #trials, propB, FR, std, relvalue]
                    newTrialTypes_all = [OVA_SO, OVB_SO, ord_SO, CJ_SO];
                    uni_newTriTypes = unique(newTrialTypes_all,'rows');
                    nnewTriTypes = size(uni_newTriTypes,1);
                    TrialTypes_tbl = uni_newTriTypes;
                    for inewTriType = 1:nnewTriTypes
                        newTriType = uni_newTriTypes(inewTriType,:);
                        ind_newTriType = ismember(newTrialTypes_all, newTriType,'rows');
                        %
                        trialnum = sum(ind_newTriType);
                        TrialTypes_tbl(inewTriType,5) = trialnum;
                        %
                        trialnum_chB = sum(CJ_SO(ind_newTriType)==-1);
                        TrialTypes_tbl(inewTriType,6) = trialnum_chB./trialnum;
                        %
                        FR_all = act_SO(ind_newTriType,:);
                        FR_itritype = nanmean(FR_all);
                        FR_std_itritype = nanstd(FR_all);
                        TrialTypes_tbl(inewTriType,7) = FR_itritype; 
                        TrialTypes_tbl(inewTriType,8) = FR_std_itritype;  
                        TrialTypes_tbl(inewTriType,9) = rho_inSO(icell);      
                    end
                    ind_badtrialtype = TrialTypes_tbl(:,5)<atleast_nntrials;
                    TrialTypes_tbl(ind_badtrialtype,:) = [];
                    %
                    % define independent variables (use same name as cell type)  
                    ABBA = TrialTypes_tbl(:,3); % -1 BA; 1 AB  
                    valA = TrialTypes_tbl(:,1);
                    valAinBA = valA;
                    % valAinBA(ABBA== 1) = nan;
                    valAinBA(ABBA== 1) = 0;
                    valAinAB = valA;
                    % valAinAB(ABBA==-1) = nan; 
                    valAinAB(ABBA==-1) = 0;
                    valB = TrialTypes_tbl(:,2);
                    valBinBA = valB;
                    % valBinBA(ABBA== 1) = nan;
                    valBinBA(ABBA== 1) = 0;
                    valBinAB = valB;
                    % valBinAB(ABBA==-1) = nan;      
                    valBinAB(ABBA==-1) = 0;   
                    val1 = valA;
                    val1(ABBA ==-1) = valB(ABBA ==-1);
                    val2 = valA;
                    val2(ABBA == 1) = valB(ABBA == 1);
                    chJ = TrialTypes_tbl(:,4); % -1 B; 1 A  
                    chV = valA;
                    chV(chJ ==-1) = valB(chJ ==-1);
                    chVA = chV;
                    % chVA(chJ ==-1) = nan;
                    chVA(chJ ==-1) = 0;
                    chVB = chV;
                    % chVB(chJ == 1) = nan;
                    chVB(chJ == 1) = 0;
                    chO = chJ.*ABBA; % -1 chosen second; 1 chosen first
                    chV1 = chV;
%                     chV1(chO ==-1) = nan;
                    chV1(chO ==-1) = 0;
                    chV2 = chV;
%                     chV2(chO == 1) = nan;  
                    chV2(chO == 1) = 0;
                    %
                    if iTW == 1
                        OVA = valAinAB;
                        OVB = valBinBA;
                        CJ = ABBA;   % use ABBA instead
                        CV = val1;   % use OV1 and OV2 instead
                        % CV = chV1;   % OV1 when 1 is chosen
                        QVA =OVA./rho_inSO(icell);
                        QVB = OVB;
                        %
                        if isequal(classname,'OVA')
                            ind_goodFR = ord_SO ==  1;  % AB trial
                        elseif isequal(classname,'OVB')
                            ind_goodFR = ord_SO == -1;  % BA trial 
                        elseif isequal(classname,'CV')
                            ind_goodFR = CO_SO == 1;  % offer1 chosen
                        else
                            ind_goodFR = logical(ones(size(act_SO))); % all trials
                        end
                    elseif iTW == 2 
                        OVA = valAinBA;
                        OVB = valBinAB;
                        CJ = ABBA;   % use ABBA instead
                        CV = val2;   % use OV1 and OV2 instead 
                        % CV = chV2;   % OV2 when 2 is chosen
                        QVA =OVA./rho_inSO(icell);
                        QVB = OVB;
                        %
                        if isequal(classname,'OVA')
                            ind_goodFR = ord_SO == -1;  % BA trial
                        elseif isequal(classname,'OVB')
                            ind_goodFR = ord_SO ==  1;  % AB trial 
                        elseif isequal(classname,'CV')
                            ind_goodFR = CO_SO == -1;  % offer2 chosen
                        else
                            ind_goodFR = logical(ones(size(act_SO))); % all trials
                        end                    
                    end
                    % fitting
                    eval(['XXX = ',classname,';'])                    
                    YYY = TrialTypes_tbl(:,7);
                    ind_good = ~isnan(XXX);
                    XXX = XXX(ind_good);
                    YYY = YYY(ind_good);
                    % XXX(isnan(XXX)) = 0;
                    g = fittype('m*x+q','coeff',{'m','q'});
                    [linfit,goodness] = fit(XXX,YYY,g);
                    %
                    FR_inSO(icell,iTW) = nanmean(act_SO(ind_goodFR));              
                    Rsq_inSO(icell,iTW) = goodness.rsquare;
                    slope_inSO(icell,iTW) = linfit.m;
                    intcept_inSO(icell,iTW) = linfit.q;
                    interval = confint(linfit,0.90);
                    if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                    nonzero_inSO(icell,iTW) = nonzero;                
                    maxintc_inSO(icell,iTW) = slope_inSO(icell,iTW).*valrange_inSO(icell,iTW)+intcept_inSO(icell,iTW);
                    %
                    % % offer value cells, QVA and QVB fitting
                    if ismember(classname,{'OVA','OVB'})
                        if iTW == 1
                            OVA = valAinAB./rho_inSO(icell);
                            OVB = valBinBA;
                        elseif iTW == 2 
                            OVA = valAinBA./rho_inSO(icell);
                            OVB = valBinAB;
                        end
                        eval(['XXX = ',classname,';'])               
                        YYY = TrialTypes_tbl(:,7);
                        ind_good = ~isnan(XXX);
                        XXX = XXX(ind_good);
                        YYY = YYY(ind_good);
                        g = fittype('m*x+q','coeff',{'m','q'});
                        [linfit,goodness] = fit(XXX,YYY,g);
                        slope_QAB_inSO(icell,iTW) = linfit.m;
                        intcept_QAB_inSO(icell,iTW) = linfit.q;
                        interval = confint(linfit,0.90);
                        if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                        nonzero_QAB_inSO(icell,iTW) =nonzero;                        
                    else
                        slope_QAB_inSO = [];
                        nonzero_QAB_inSO = [];
                        intcept_QAB_inSO = [];
                    end
                    
                end % for iTW                              
            end %for icell
            %
            eval(['allcellnames.',classname,'.',slopesignname,' = cellnames_iclass;'])
            eval(['allmonkeys.',classname,'.',slopesignname,' = monkeyname_list;'])
            %
            eval(['allFRs.',classname,'.',slopesignname,'.JC = FR_inJC;'])
            eval(['allFRs.',classname,'.',slopesignname,'.SO = FR_inSO;'])
            eval(['allSteepness.',classname,'.',slopesignname,'.JC = steepness_inJC;'])
            eval(['allSteepness.',classname,'.',slopesignname,'.SO = steepness_inSO;'])
            eval(['allOrderbias.',classname,'.',slopesignname,'.SO = orderbias_inSO;'])
            eval(['allsessrange.',classname,'.',slopesignname,'.JC = sessrange_inJC;'])
            eval(['allsessrange.',classname,'.',slopesignname,'.SO = sessrange_inSO;'])            
            eval(['allrhos.',classname,'.',slopesignname,'.JC = rho_inJC;'])
            eval(['allrhos.',classname,'.',slopesignname,'.SO = rho_inSO;'])
            %
            eval(['allslopes.',classname,'.',slopesignname,'.JC = slope_inJC;'])
            eval(['allslopes.',classname,'.',slopesignname,'.SO = slope_inSO;'])           
            eval(['allRsq.',classname,'.',slopesignname,'.JC = Rsq_inJC;'])
            eval(['allRsq.',classname,'.',slopesignname,'.SO = Rsq_inSO;'])           
            eval(['allnonzeros.',classname,'.',slopesignname,'.JC = nonzero_inJC;'])
            eval(['allnonzeros.',classname,'.',slopesignname,'.SO = nonzero_inSO;'])
            eval(['allintcepts.',classname,'.',slopesignname,'.JC = intcept_inJC;'])
            eval(['allintcepts.',classname,'.',slopesignname,'.SO = intcept_inSO;'])
            eval(['allmaxinters.',classname,'.',slopesignname,'.JC = maxintc_inJC;'])
            eval(['allmaxinters.',classname,'.',slopesignname,'.SO = maxintc_inSO;'])            
            %
            eval(['allslopes_QAB.',classname,'.',slopesignname,'.JC = slope_QAB_inJC;'])            
            eval(['allintcepts_QAB.',classname,'.',slopesignname,'.JC = intcept_QAB_inJC;'])    
            eval(['allnonzeros_QAB.',classname,'.',slopesignname,'.JC = nonzero_QAB_inJC;']) 
            eval(['allslopes_QAB.',classname,'.',slopesignname,'.SO = slope_QAB_inSO;'])            
            eval(['allintcepts_QAB.',classname,'.',slopesignname,'.SO = intcept_QAB_inSO;'])    
            eval(['allnonzeros_QAB.',classname,'.',slopesignname,'.SO = nonzero_QAB_inSO;']) 
           
        end % try catch end
        
        if ~slopepooled
            plot_pop_preference_bias_JCSOinTT
        end
        
    end %for islopesign
end %for iclassname

% %
% filename = ['pop_encoding_slope_JCSOinTT_redoslope',differentRho,'_',monkey_ana];
filename = ['pop_preference_bias_JCSOinTT_redoslope',differentRho,'_',monkey_ana, dobhvlogit]; 
eval(['save ',filename ' allcellnames allmonkeys allFRs allSteepness allOrderbias allsessrange allrhos '...
                       ' allslopes allnonzeros allRsq allintcepts allmaxinters '...
                       ' allslopes_QAB  allintcepts_QAB allnonzeros_QAB '...
                       ])

% %
% pool positive and negative : only for OV cells
% %
if slopepooled
    cellclasses = {'OVA','OVB'};
    ncellclasses = size(cellclasses,2);
    n = size(slopesigns,2);
    for icellclass = 1:ncellclasses
        classname = cellclasses{icellclass};                    
        slopesignname = 'slopemerged';
        %
        eval(['cellnames_iclass = {allcellnames.',classname,'.positive{:}, allcellnames.',classname,'.negative{:}}'';'])
        eval(['monkeyname_list = {allmonkeys.',classname,'.positive{:}, allmonkeys.',classname,'.negative{:}}'';'])
        %
        eval(['steepness_inJC = [allSteepness.',classname,'.positive.JC; allSteepness.',classname,'.negative.JC];'])
        eval(['steepness_inSO = [allSteepness.',classname,'.positive.SO; allSteepness.',classname,'.negative.SO];'])
        eval(['orderbias_inSO = [allOrderbias.',classname,'.positive.SO; allOrderbias.',classname,'.negative.SO];'])
        eval(['sessrange_inJC = [allsessrange.',classname,'.positive.JC; allsessrange.',classname,'.negative.JC];'])
        eval(['sessrange_inSO = [allsessrange.',classname,'.positive.SO; allsessrange.',classname,'.negative.SO];'])
        eval(['rho_inJC = [allrhos.',classname,'.positive.JC; allrhos.',classname,'.negative.JC];'])
        eval(['rho_inSO = [allrhos.',classname,'.positive.SO; allrhos.',classname,'.negative.SO];'])
        %
        eval(['FR_inJC = [allFRs.',classname,'.positive.JC; allFRs.',classname,'.negative.JC];'])
        eval(['FR_inSO = [allFRs.',classname,'.positive.SO; allFRs.',classname,'.negative.SO];'])
        eval(['nonzero_inJC = [allnonzeros.',classname,'.positive.JC; allnonzeros.',classname,'.negative.JC];'])
        eval(['nonzero_inSO = [allnonzeros.',classname,'.positive.SO; allnonzeros.',classname,'.negative.SO];'])
        eval(['slope_inJC = [allslopes.',classname,'.positive.JC; -allslopes.',classname,'.negative.JC];'])
        eval(['slope_inSO = [allslopes.',classname,'.positive.SO; -allslopes.',classname,'.negative.SO];'])
        eval(['intcept_inJC = [allintcepts.',classname,'.positive.JC; -allintcepts.',classname,'.negative.JC];'])
        eval(['intcept_inSO = [allintcepts.',classname,'.positive.SO; -allintcepts.',classname,'.negative.SO];'])
        eval(['maxintc_inJC = [allmaxinters.',classname,'.positive.JC; allmaxinters.',classname,'.negative.JC];'])
        eval(['maxintc_inSO = [allmaxinters.',classname,'.positive.SO; allmaxinters.',classname,'.negative.SO];'])
        eval(['Rsq_inJC = [allRsq.',classname,'.positive.JC; allRsq.',classname,'.negative.JC];'])
        eval(['Rsq_inSO = [allRsq.',classname,'.positive.SO; allRsq.',classname,'.negative.SO];'])   
        %
        eval(['nonzero_QAB_inJC = [allnonzeros_QAB.',classname,'.positive.JC; allnonzeros_QAB.',classname,'.negative.JC];'])
        eval(['nonzero_QAB_inSO = [allnonzeros_QAB.',classname,'.positive.SO; allnonzeros_QAB.',classname,'.negative.SO];'])
        eval(['slope_QAB_inJC = [allslopes_QAB.',classname,'.positive.JC; -allslopes_QAB.',classname,'.negative.JC];'])
        eval(['slope_QAB_inSO = [allslopes_QAB.',classname,'.positive.SO; -allslopes_QAB.',classname,'.negative.SO];'])
        eval(['intcept_QAB_inJC = [allintcepts_QAB.',classname,'.positive.JC; -allintcepts_QAB.',classname,'.negative.JC];'])
        eval(['intcept_QAB_inSO = [allintcepts_QAB.',classname,'.positive.SO; -allintcepts_QAB.',classname,'.negative.SO];'])
        
        actrange_inJC = abs(maxintc_inJC - intcept_inJC);
        actrange_inSO = abs(maxintc_inSO - intcept_inSO);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % %
        ind_G = ismember(monkeyname_list,'G');
        ind_J = ismember(monkeyname_list,'J');
        % remove outlier based on steepness; kout IQR methods
        kout = 1.5;
        % % Gerinvho
%         steepness_inJC_G =  steepness_inJC(ind_G);
%         quantiles_JC = quantile(steepness_inJC_G,[0.25 0.5 0.75]);
%         IQR_JC = quantiles_JC(3) - quantiles_JC(1);
%         steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
%         steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
%         %
%         steepness_inSO_G =  steepness_inSO(ind_G);
%         quantiles_SO = quantile(steepness_inSO_G,[0.25 0.5 0.75]);
%         IQR_SO = quantiles_SO(3) - quantiles_SO(1);
%         steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];  
%         steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
%         %
%         deltasteep = steepness_inJC - steepness_inSO;
%         deltasteep_G = deltasteep(ind_G);
%         quantiles_del = quantile(deltasteep_G,[0.25 0.5 0.75]);
%         IQR_del = quantiles_del(3) - quantiles_del(1);
%         steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
%         steepoutlier_del(2) = min([steepoutlier_del(2),20]);
%         %
%         ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_G;
%         ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_G;
%         ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_G;
%         ind_noneout_G = ind_goodJC & ind_goodSO; % & ind_gooddelta;
        ind_goodJC = (steepness_inJC > steepoutlier_JC_G(1)  & steepness_inJC < steepoutlier_JC_G(2)) & ind_G;
        ind_goodSO = (steepness_inSO > steepoutlier_SO_G(1)  & steepness_inSO < steepoutlier_SO_G(2)) & ind_G;
        ind_noneout_G = ind_goodJC & ind_goodSO; % & ind_gooddelta;
        % % Juan
%         steepness_inJC_J =  steepness_inJC(ind_J);
%         quantiles_JC = quantile(steepness_inJC_J,[0.25 0.5 0.75]);
%         IQR_JC = quantiles_JC(3) - quantiles_JC(1);
%         steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
%         steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
%         %
%         steepness_inSO_J =  steepness_inSO(ind_J);
%         quantiles_SO = quantile(steepness_inSO_J,[0.25 0.5 0.75]);
%         IQR_SO = quantiles_SO(3) - quantiles_SO(1);
%         steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];   
%         steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
%         %
%         deltasteep = steepness_inJC - steepness_inSO;
%         deltasteep_J = deltasteep(ind_J);
%         quantiles_del = quantile(deltasteep_J,[0.25 0.5 0.75]);
%         IQR_del = quantiles_del(3) - quantiles_del(1);
%         steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
%         steepoutlier_del(2) = min([steepoutlier_del(2),20]);
%         %
%         ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_J;
%         ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_J;
%         ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_J;
        ind_goodJC = (steepness_inJC > steepoutlier_JC_J(1)  & steepness_inJC < steepoutlier_JC_J(2)) & ind_J;
        ind_goodSO = (steepness_inSO > steepoutlier_SO_J(1)  & steepness_inSO < steepoutlier_SO_J(2)) & ind_J;
        ind_noneout_J = ind_goodJC & ind_goodSO; % & ind_gooddelta;
        % %
        ind_noneout = ind_noneout_G | ind_noneout_J;
        % %
        % ind_rightslope = slope_inJC>0 & slope_inSO>0;
        % %
        if ~dosteepout
            ind_noneout = logical(ones(size(ind_noneout)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
       
        %
        plot_pop_preference_bias_JCSOinTT
    end
end




% %
% save figures
% %
if saveplots
h =  findobj('type','figure');
nplots = length(h);
for iplot = 1:nplots
    hf = figure(iplot);
    set(gcf,'Units','normalized', 'position',[0.5 0 0.5 1], 'PaperPositionMode','auto')
    eval(['print -fillpage -opengl -dpsc2 -append -r0 ',figuresave]);
    close
end
end



%%
function [Rsqall, Rsq_max, ind_max] = find_Rsq(cellname,cellclass, JCorSO, JCSOclassifySep)
session = cellname(1:8);
readsession_TT
filename = [dirroot,cellname,'_cellstats'];
eval(['load ',filename])
%

cellclass_all = {'OVA','OVB','CV','CJ'};
ind_class = find(ismember(cellclass_all,cellclass)==1);

if isequal(JCorSO,'JC')
    pairname = 'AB';
    if ~JCSOclassifySep
        TWall = [2, 3, 7];
        modelmatrix = [12 12 15;
                       13 13 16;
                        6  6  6;
                       14 14 14; ];
    elseif JCSOclassifySep
        TWall = [2, 3, 6, 7];
        modelmatrix = [12 12 12 15;
                       13 13 13 16;
                        6  6  6  6;
                       14 14 14 14; ];
    end       
elseif isequal(JCorSO,'SO')
    pairname = 'ABA';
    TWall = [2, 4, 8];
    modelmatrix = [ 2  3 16;
                    6  5 17;
                    7  9 15;
                   20 20 19; ];
end

modelnum = modelmatrix(ind_class,:);

eval(['Rsq_tgt = diag(cellstats.tuningfit.',JCorSO,'.',pairname,'.Rsq(modelnum,TWall));'])

[Rsq_max, ind_max] = max(Rsq_tgt);

Rsqall = Rsq_tgt;

end

%%
function plotErrorEllipse(mu_ell, Sigma_ell, p_ell)
    s = -2 * log(1 - p_ell);
    [V, D] = eig(Sigma_ell * s);
    t = linspace(0, 2 * pi);
    a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
    plot(a(1, :) + mu_ell(1), a(2, :) + mu_ell(2), 'k-');
end
