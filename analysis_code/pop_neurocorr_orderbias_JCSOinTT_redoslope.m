% pop_neurocorr_orderbias_JCSOinTT.m
%
% This script analyzes the neural correlates for order bias in Task 2 (SO).
% This script focuses on postoffer time window for JC and postoffer1 and postoffer2 time window for SO
% offer value cells  - EO (offer1 TW) v.s. OE (offer2 TW): slope and FR
% chosen value cells - offer1 TW v.s. offer2 TW: slope and FR 
%                    - rho_neuronal related comparison (in analysis: pop_rho_neuronal_JCSOinTT.m)
% chosen juice cells - circuit inhibition related comparison (in analysis: circuit_inhibition_xx_TT.m)
% 
% This script uses plot_pop_neurocorr_orderbias_JCSOinTT.m to plot 
%

%
% author:     Nov 2019; WS
% revision:   Feb 2020; WS - add Rsq 
% revision:   Aug 2020: WS - add full slope and redo the regression
% revision:   Oct 2020: WS - finalize the analysis, fix the problem
% revision:   Jan 2021: WS - add CV cells and CJ cells; add postjuice time
%                            window
% revision:   Jan 2021: WS - focus on order bias

close all
clearvars

brainarea = 'OFC'; % 'DLPFC', 'VLPFC' 'OFC'
monkey_ana = 'both'; % 'Gervinho', 'Juan', 'both'

% % % 
JCSOclassifySep = 0; % if 1, load neurons that is classified seperately by JC or SO trials: 0-JCSOcellist

juiceABpooled = 1;   % if 1, pool cell groups: A and B
slopepooled = 1;     % if 1, pool cell postive and negative

atleast_nntrials = 2;

dobhvlogit = ''; % '': probit; '_logit'; '_neworderbias'

saveplots = 0;
if saveplots
if      juiceABpooled &  slopepooled
    figuresave = ['C:\Experiments\TwoTasks\Analysis\Analysis_OFC\pop_neurocorr_orderbias_JCSOinTT_redoslope_juiceABpooled_slopepooled_',monkey_ana];
elseif  juiceABpooled & ~slopepooled
    figuresave = ['C:\Experiments\TwoTasks\Analysis\Analysis_OFC\pop_neurocorr_orderbias_JCSOinTT_redoslope_juiceABpooled_',monkey_ana];
elseif ~juiceABpooled &  slopepooled
    figuresave = ['C:\Experiments\TwoTasks\Analysis\Analysis_OFC\pop_neurocorr_orderbias_JCSOinTT_redoslope_slopepooled_',monkey_ana];
elseif ~juiceABpooled & ~slopepooled
    figuresave = ['C:\Experiments\TwoTasks\Analysis\Analysis_OFC\pop_neurocorr_orderbias_JCSOinTT_redoslope_',monkey_ana];    
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
allslopes     = [];
allnonzeros   = [];
allRsq        = [];
allintcepts   = [];
allmaxinters  = [];
allFRs        = [];
allSteepness  = [];
allOrderbias  = [];
allsessrange  = [];
allmonkeys    = [];
allslopes_QAB = [];     % for CV cells
allslopes_chV12 = [];   % for CV cells, postjuice TW 
allintcepts_chV12 = []; % for CV cells, postjuice TW 
allnonzeros_QAB = [];   % for CV cells
allnonzeros_chV12 = []; % for CV cells, postjuice TW  
allrhos       = [];

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
% cellclassnames = {'OVA', 'OVB', 'CV', 'CJ'};
cellclassnames = {'OVA', 'OVB', 'CV'};
% cellclassnames = {'OVA', 'OVB'};
% cellclassnames = {'CV'};
% cellclassnums = [1,2,3,4];
cellclassnums = [1,2, 3];
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
            % filename = ['pop_encoding_slope_JCSOinTT_redoslope',differentRho,'_',monkey_ana];
            filename = ['pop_neurocorr_orderbias_JCSOinTT_redoslope',differentRho,'_',monkey_ana,dobhvlogit];            
            load(filename);            
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
            eval(['FR_inJC = allFRs.',classname,'.',slopesignname,'.JC;'])
            eval(['FR_inSO = allFRs.',classname,'.',slopesignname,'.SO;'])
            eval(['steepness_inJC = allSteepness.',classname,'.',slopesignname,'.JC;'])
            eval(['steepness_inSO = allSteepness.',classname,'.',slopesignname,'.SO;'])
            eval(['orderbias_inSO = allOrderbias.',classname,'.',slopesignname,'.SO;'])
            eval(['sessrange_inJC = allsessrange.',classname,'.',slopesignname,'.JC;'])
            eval(['sessrange_inSO = allsessrange.',classname,'.',slopesignname,'.SO;'])
            eval(['cellnames_iclass = allcellnames.',classname,'.',slopesignname,';'])
            eval(['monkeyname_list = allmonkeys.',classname,'.',slopesignname,';'])
            eval(['rho_inJC = allrhos.',classname,'.',slopesignname,'.JC;'])
            eval(['rho_inSO = allrhos.',classname,'.',slopesignname,'.SO;'])
            if isequal(classname,'CV')
                eval(['slopes_QAB_inSO = allslopes_QAB.',classname,'.',slopesignname,'.SO;'])
                eval(['slopes_chV12_inSO = allslopes_chV12.',classname,'.',slopesignname,'.SO;'])
                eval(['intcepts_chV12_inSO = allintcepts_chV12.',classname,'.',slopesignname,'.SO;'])
                eval(['nonzeros_QAB_inSO = allnonzeros_QAB.',classname,'.',slopesignname,'.SO;']) 
                eval(['nonzeros_chV12_inSO = allnonzeros_chV12.',classname,'.',slopesignname,'.SO;']) 
            end
            
        catch
            
            ind_JCiclass = JCclasses==cellclassnums(iclassname) & JCslopes==slopesigns(islopesign);
            ind_SOiclass = SOclasses==cellclassnums(iclassname) & SOslopes==slopesigns(islopesign);
            ind_iclass = ind_JCiclass & ind_SOiclass;
            % ind_iclass = ind_SOiclass;
            %
            ncells = sum(ind_iclass);
            cellnames_iclass = cellnames(ind_iclass);
            cellnum_iclass = find(ind_iclass==1);

            monkeyname_list = [];

            slope_inJC = [];     slope_inSO = [];
            Rsq_inJC = [];       Rsq_inSO = [];
            nonzero_inJC = [];   nonzero_inSO = [];
            FR_inJC = [];        FR_inSO = [];
            intcept_inJC = [];   intcept_inSO = [];
            maxintc_inJC = [];   maxintc_inSO = [];
            valrange_inJC = [];  valrange_inSO = [];
            steepness_inJC = []; steepness_inSO = [];
            rho_inJC = [];       rho_inSO = [];
            sessrange_inJC = []; sessrange_inSO = [];
            orderbias_inSO = [];
            
            slopes_QAB_inSO = [];
            slopes_chV12_inSO = [];
            intcepts_chV12_inSO = [];
            nonzeros_QAB_inSO = [];
            nonzeros_chV12_inSO = [];
            
            
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
                        % orderbias_inSO(icell,:) = -2.*rho_inSO(icell,:).*psyphycell.SO.OrdChHyst.sigmoidfit.beta(3)./psyphycell.SO.OrdChHyst.sigmoidfit.beta(2);
                    end
                     
                end                
                %
                
                % % % 
                filename = [dirroot,cellname,'_tuning'];
                eval(['load ',filename])
                
                % redo the linear regression
                % three time windows: JC: postoffer  postoffer  postjuice               
                twins = {'postoffer','postoffer','postjuice'};
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
                    % valrange_CV_JC = nanmax(CV_JC) - nanmin(CV_JC);
                    valrange_CV_JC = nanmax(CV_JC);
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
                    % FR_inJC(icell,iTW) = nanmean(act_JC);    
                    FR_inJC(icell,iTW) = linfit.q + linfit.m *valrange_inJC(icell,iTW)/2;
                    Rsq_inJC(icell,iTW) = goodness.rsquare;
                    slope_inJC(icell,iTW) = linfit.m;
                    intcept_inJC(icell,iTW) = linfit.q;
                    interval = confint(linfit,0.90);
                    if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                    nonzero_inJC(icell,iTW) = nonzero;                
                    maxintc_inJC(icell,iTW) = slope_inJC(icell,iTW).*valrange_inJC(icell,iTW)+intcept_inJC(icell,iTW);
                end
                
                % %
                % redo the linear regression
                % three time windows: SO: postoffer1 postoffer2 postjuice
                twins = {'postoffer1','postoffer2','postjuice'};
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
                    % valrange_CV_SO = nanmax(CV_SO) - nanmin(CV_SO);
                    valrange_CV_SO = nanmax(CV_SO);
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
%                     valAinBA(ABBA== 1) = nan;
                    valAinBA(ABBA== 1) = 0;
                    valAinAB = valA;
%                     valAinAB(ABBA==-1) = nan; 
                    valAinAB(ABBA==-1) = 0;
                    valB = TrialTypes_tbl(:,2);
                    valBinBA = valB;
%                     valBinBA(ABBA== 1) = nan;
                    valBinBA(ABBA== 1) = 0;
                    valBinAB = valB;
%                     valBinAB(ABBA==-1) = nan;      
                    valBinAB(ABBA==-1) = 0;   
                    val1 = valA;
                    val1(ABBA ==-1) = valB(ABBA ==-1);
                    val2 = valA;
                    val2(ABBA == 1) = valB(ABBA == 1);
                    chJ = TrialTypes_tbl(:,4); % -1 B; 1 A  
                    chV = valA;
                    chV(chJ ==-1) = valB(chJ ==-1);
                    chVA = chV;
%                     chVA(chJ ==-1) = nan;
                    chVA(chJ ==-1) = 0;
                    chVB = chV;
%                     chVB(chJ == 1) = nan;
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
                    elseif iTW == 3 
                        OVA = chVA;
                        OVB = chVB;
                        CJ = chJ;   
                        CV = chV;  
                        QVA =OVA./rho_inSO(icell);
                        QVB = OVB;  
                        %
                        if isequal(classname,'OVA')
                            ind_goodFR = CJ_SO ==  1;  % A chosen
                        elseif isequal(classname,'OVB')
                            ind_goodFR = CJ_SO == -1;  % B chosen 
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
                    % FR_inSO(icell,iTW) = nanmean(act_SO(ind_goodFR));
                    FR_inSO(icell,iTW) = linfit.q + linfit.m *valrange_inSO(icell,iTW)/2;
                    Rsq_inSO(icell,iTW) = goodness.rsquare;
                    slope_inSO(icell,iTW) = linfit.m;
                    intcept_inSO(icell,iTW) = linfit.q;
                    interval = confint(linfit,0.90);
                    if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                    nonzero_inSO(icell,iTW) = nonzero;                
                    maxintc_inSO(icell,iTW) = slope_inSO(icell,iTW).*valrange_inSO(icell,iTW)+intcept_inSO(icell,iTW);
                    %
                    % % CV cells, QVA and QVB fitting
                    if isequal(classname,'CV')
                        XXX = QVA;                   
                        YYY = TrialTypes_tbl(:,7);
                        ind_good = ~isnan(XXX);
                        XXX = XXX(ind_good);
                        YYY = YYY(ind_good);
                        g = fittype('m*x+q','coeff',{'m','q'});
                        [linfit,goodness] = fit(XXX,YYY,g);
                        slopes_QAB_inSO(icell,iTW) = linfit.m;
                        interval = confint(linfit,0.90);
                        if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                        nonzeros_QAB_inSO(icell,iTW) =nonzero;
                        %
                        XXX = QVB;                   
                        YYY = TrialTypes_tbl(:,7);
                        ind_good = ~isnan(XXX);
                        XXX = XXX(ind_good);
                        YYY = YYY(ind_good);
                        g = fittype('m*x+q','coeff',{'m','q'});
                        [linfit,goodness] = fit(XXX,YYY,g);
                        slopes_QAB_inSO(icell,iTW+3) = linfit.m;
                        interval = confint(linfit,0.90);
                        if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                        nonzeros_QAB_inSO(icell,iTW+3) =  nonzero;
                    else
                        slopes_QAB_inSO = [];
                        nonzeros_QAB_inSO = [];
                    end
                    %
                    % % CV cells, postjuice TW, chosen value 1 and chosen value 2
                    if isequal(classname,'CV')
                    if iTW == 3 
                        XXX = chV1;                   
                        YYY = TrialTypes_tbl(:,7);
                        ind_good = ~isnan(XXX);
                        XXX = XXX(ind_good);
                        YYY = YYY(ind_good);
                        g = fittype('m*x+q','coeff',{'m','q'});
                        [linfit,goodness] = fit(XXX,YYY,g);
                        slopes_chV12_inSO(icell,1) = linfit.m;
                        intcepts_chV12_inSO(icell,1) = linfit.q;
                        interval = confint(linfit,0.90);
                        if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                        nonzeros_chV12_inSO(icell,1) = nonzero;
                        %
                        XXX = chV2;                   
                        YYY = TrialTypes_tbl(:,7);
                        ind_good = ~isnan(XXX);
                        XXX = XXX(ind_good);
                        YYY = YYY(ind_good);
                        g = fittype('m*x+q','coeff',{'m','q'});
                        [linfit,goodness] = fit(XXX,YYY,g);
                        slopes_chV12_inSO(icell,2) = linfit.m;
                        intcepts_chV12_inSO(icell,2) = linfit.q;
                        interval = confint(linfit,0.90);
                        if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                        nonzeros_chV12_inSO(icell,2) = nonzero;
                    end
                    else
                        slopes_chV12_inSO = [];
                        intcepts_chV12_inSO = [];
                        nonzeros_chV12_inSO = [];
                    end
                    
                end % for iTW                              
            end %for icell
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
            eval(['allFRs.',classname,'.',slopesignname,'.JC = FR_inJC;'])
            eval(['allFRs.',classname,'.',slopesignname,'.SO = FR_inSO;'])
            eval(['allSteepness.',classname,'.',slopesignname,'.JC = steepness_inJC;'])
            eval(['allSteepness.',classname,'.',slopesignname,'.SO = steepness_inSO;'])
            eval(['allOrderbias.',classname,'.',slopesignname,'.SO = orderbias_inSO;'])
            eval(['allsessrange.',classname,'.',slopesignname,'.JC = sessrange_inJC;'])
            eval(['allsessrange.',classname,'.',slopesignname,'.SO = sessrange_inSO;'])
            eval(['allcellnames.',classname,'.',slopesignname,' = cellnames_iclass;'])
            eval(['allmonkeys.',classname,'.',slopesignname,' = monkeyname_list;'])
            eval(['allrhos.',classname,'.',slopesignname,'.JC = rho_inJC;'])
            eval(['allrhos.',classname,'.',slopesignname,'.SO = rho_inSO;'])
            eval(['allslopes_QAB.',classname,'.',slopesignname,'.SO = slopes_QAB_inSO;']) 
            eval(['allslopes_chV12.',classname,'.',slopesignname,'.SO = slopes_chV12_inSO;']) 
            eval(['allintcepts_chV12.',classname,'.',slopesignname,'.SO = intcepts_chV12_inSO;']) 
   
            eval(['allnonzeros_QAB.',classname,'.',slopesignname,'.SO = nonzeros_QAB_inSO;']) 
            eval(['allnonzeros_chV12.',classname,'.',slopesignname,'.SO = nonzeros_chV12_inSO;']) 
            
            
            
        end % try catch end
        
        % % % 
        % plot each cell class
        % % %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        if ~juiceABpooled & ~slopepooled
        %
        ind_G = ismember(monkeyname_list,'G')';
        ind_J = ismember(monkeyname_list,'J')';
        % remove outlier based on steepness; kout IQR methods
        kout = 1.5;
        % % Gerinvho
        steepness_inJC_G =  steepness_inJC(ind_G);
        quantiles_JC = quantile(steepness_inJC_G,[0.25 0.5 0.75]);
        IQR_JC = quantiles_JC(3) - quantiles_JC(1);
        steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
        steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
        %
        steepness_inSO_G =  steepness_inSO(ind_G);
        quantiles_SO = quantile(steepness_inSO_G,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];     
        steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
        %
        orderbias_inSO_G =  orderbias_inSO(ind_G);
        quantiles_SO = quantile(orderbias_inSO_G,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        orderoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];             
        %
        deltasteep = steepness_inJC - steepness_inSO;
        deltasteep_G = deltasteep(ind_G);
        quantiles_del = quantile(deltasteep_G,[0.25 0.5 0.75]);
        IQR_del = quantiles_del(3) - quantiles_del(1);
        steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
        steepoutlier_del(2) = min([steepoutlier_del(2),20]);
        %
        ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_G;
        ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_G;
        ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_G;
        ind_noneout_G = ind_goodJC & ind_goodSO; % & ind_gooddelta;       
        ind_goodorder_G = (orderbias_inSO > orderoutlier_SO(1)  & orderbias_inSO < orderoutlier_SO(2)) & ind_G;
        % % 
        % % Juan
        steepness_inJC_J =  steepness_inJC(ind_J);
        quantiles_JC = quantile(steepness_inJC_J,[0.25 0.5 0.75]);
        IQR_JC = quantiles_JC(3) - quantiles_JC(1);
        steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
        steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
        %
        steepness_inSO_J =  steepness_inSO(ind_J);
        quantiles_SO = quantile(steepness_inSO_J,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];  
        steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
        %
        orderbias_inSO_J =  orderbias_inSO(ind_J);
        quantiles_SO = quantile(orderbias_inSO_J,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        orderoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];
        %
        deltasteep = steepness_inJC - steepness_inSO;
        deltasteep_J = deltasteep(ind_J);
        quantiles_del = quantile(deltasteep_J,[0.25 0.5 0.75]);
        IQR_del = quantiles_del(3) - quantiles_del(1);
        steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
        steepoutlier_del(2) = min([steepoutlier_del(2),20]);
        %
        ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_J;
        ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_J;
        ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_J;
        ind_noneout_J = ind_goodJC & ind_goodSO; % & ind_gooddelta;  
        ind_goodorder_J = (orderbias_inSO > orderoutlier_SO(1)  & orderbias_inSO < orderoutlier_SO(2)) & ind_J;
        % %
        ind_noneout = ind_noneout_G | ind_noneout_J;
        ind_goodorder = ind_goodorder_G | ind_goodorder_J;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
         
%         % % %
%         % ANOVA threshold
%         ind_ANOVA_JC  = pvalANOVA_inJC(:,2)<0.01;
%         ind_ANOVA_SO1 = pvalANOVA_inSO(:,2)<0.01;
%         ind_ANOVA_SO2 = pvalANOVA_inSO(:,4)<0.01;
        
        %
        actrange_inJC = abs(maxintc_inJC - intcept_inJC);
        actrange_inSO = abs(maxintc_inSO - intcept_inSO);
        %
        plot_pop_neurocorr_orderbias_JCSOinTT
        
        end
    end %for islopesign
end %for iclassname

% %
% filename = ['pop_encoding_slope_JCSOinTT_redoslope',differentRho,'_',monkey_ana];
filename = ['pop_neurocorr_orderbias_JCSOinTT_redoslope',differentRho,'_',monkey_ana,dobhvlogit];
eval(['save ',filename ' allrhos allslopes allnonzeros allRsq allintcepts allmaxinters '...
                       ' allFRs allSteepness allOrderbias allsessrange allcellnames allmonkeys '...
                       ' allslopes_QAB allslopes_chV12 allintcepts_chV12 allnonzeros_QAB allnonzeros_chV12 '])


% %
% pool cell groups A and B
% %
if juiceABpooled & ~slopepooled
%     cellclasses = {'OV','CV','CJ'};
    cellclasses = {'OV'};
    ncellclasses = size(cellclasses,2);
    slopesignnames = {'positive','negative'};
    nslopesigns = size(slopesigns,2);
    for icellclass = 1:ncellclasses
        classname = cellclasses{icellclass};             
        
        for islopesign = 1:nslopesigns
            slopesignname = slopesignnames{islopesign};
        
            if isequal(classname,'OV')
                eval(['FR_inJC = [allFRs.OVA.',slopesignname,'.JC; allFRs.OVB.',slopesignname,'.JC];'])
                eval(['FR_inSO = [allFRs.OVA.',slopesignname,'.SO; allFRs.OVB.',slopesignname,'.SO];'])
                eval(['slope_inJC = [allslopes.OVA.',slopesignname,'.JC; allslopes.OVB.',slopesignname,'.JC];'])
                eval(['slope_inSO = [allslopes.OVA.',slopesignname,'.SO; allslopes.OVB.',slopesignname,'.SO];'])             
                eval(['nonzero_inJC = [allnonzeros.OVA.',slopesignname,'.JC; allnonzeros.OVB.',slopesignname,'.JC];'])
                eval(['nonzero_inSO = [allnonzeros.OVA.',slopesignname,'.SO; allnonzeros.OVB.',slopesignname,'.SO];'])
                eval(['intcept_inJC = [allintcepts.OVA.',slopesignname,'.JC; allintcepts.OVB.',slopesignname,'.JC];'])
                eval(['intcept_inSO = [allintcepts.OVA.',slopesignname,'.SO; allintcepts.OVB.',slopesignname,'.SO];'])
                eval(['maxintc_inJC = [allmaxinters.OVA.',slopesignname,'.JC; allmaxinters.OVB.',slopesignname,'.JC];'])
                eval(['maxintc_inSO = [allmaxinters.OVA.',slopesignname,'.SO; allmaxinters.OVB.',slopesignname,'.SO];'])
                eval(['steepness_inJC = [allSteepness.OVA.',slopesignname,'.JC; allSteepness.OVB.',slopesignname,'.JC];'])
                eval(['steepness_inSO = [allSteepness.OVA.',slopesignname,'.SO; allSteepness.OVB.',slopesignname,'.SO];'])
                eval(['orderbias_inSO = [allOrderbias.OVA.',slopesignname,'.SO; allOrderbias.OVB.',slopesignname,'.SO];'])
                eval(['sessrange_inJC = [allsessrange.OVA.',slopesignname,'.JC; allsessrange.OVB.',slopesignname,'.JC];'])
                eval(['sessrange_inSO = [allsessrange.OVA.',slopesignname,'.SO; allsessrange.OVB.',slopesignname,'.SO];'])
                eval(['Rsq_inJC = [allRsq.OVA.',slopesignname,'.JC; allRsq.OVB.',slopesignname,'.JC];'])
                eval(['Rsq_inSO = [allRsq.OVA.',slopesignname,'.SO; allRsq.OVB.',slopesignname,'.SO];'])           
                eval(['cellnames_iclass = {allcellnames.OVA.',slopesignname,'{:}, allcellnames.OVB.',slopesignname,'{:}}'';'])
                eval(['monkeyname_list = {allmonkeys.OVA.',slopesignname,'{:}, allmonkeys.OVB.',slopesignname,'{:}}'';'])           
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % %
            ind_G = ismember(monkeyname_list,'G');
            ind_J = ismember(monkeyname_list,'J');
            % remove outlier based on steepness; kout IQR methods
            kout = 1.5;
            % % Gerinvho
            steepness_inJC_G =  steepness_inJC(ind_G);
            quantiles_JC = quantile(steepness_inJC_G,[0.25 0.5 0.75]);
            IQR_JC = quantiles_JC(3) - quantiles_JC(1);
            steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
            steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
            %
            steepness_inSO_G =  steepness_inSO(ind_G);
            quantiles_SO = quantile(steepness_inSO_G,[0.25 0.5 0.75]);
            IQR_SO = quantiles_SO(3) - quantiles_SO(1);
            steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];  
            steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
            %
            deltasteep = steepness_inJC - steepness_inSO;
            deltasteep_G = deltasteep(ind_G);
            quantiles_del = quantile(deltasteep_G,[0.25 0.5 0.75]);
            IQR_del = quantiles_del(3) - quantiles_del(1);
            steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
            steepoutlier_del(2) = min([steepoutlier_del(2),20]);
            %
            ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_G;
            ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_G;
            ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_G;
            ind_noneout_G = ind_goodJC & ind_goodSO; % & ind_gooddelta;
            % % Juan
            steepness_inJC_J =  steepness_inJC(ind_J);
            quantiles_JC = quantile(steepness_inJC_J,[0.25 0.5 0.75]);
            IQR_JC = quantiles_JC(3) - quantiles_JC(1);
            steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
            steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
            %
            steepness_inSO_J =  steepness_inSO(ind_J);
            quantiles_SO = quantile(steepness_inSO_J,[0.25 0.5 0.75]);
            IQR_SO = quantiles_SO(3) - quantiles_SO(1);
            steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];  
            steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
            %
            deltasteep = steepness_inJC - steepness_inSO;
            deltasteep_J = deltasteep(ind_J);
            quantiles_del = quantile(deltasteep_J,[0.25 0.5 0.75]);
            IQR_del = quantiles_del(3) - quantiles_del(1);
            steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
            steepoutlier_del(2) = min([steepoutlier_del(2),20]);
            %
            ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_J;
            ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_J;
            ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_J;
            ind_noneout_J = ind_goodJC & ind_goodSO; % & ind_gooddelta;
            % %
            ind_noneout = ind_noneout_G | ind_noneout_J;
            %%%%%%%%%%%%%%%%%%%%%%%%%
           
            %
            plot_pop_neurocorr_orderbias_JCSOinTT
            
        end % for islopesign
    end % for icellclass
end % if


% %
% pool positive and negative : only for OV cells
% %
if slopepooled & ~juiceABpooled
    cellclasses = {'OVA','OVB'};
    ncellclasses = size(cellclasses,2);
    n = size(slopesigns,2);
    for icellclass = 1:ncellclasses
        classname = cellclasses{icellclass};                    
        slopesignname = 'slopemerged';
        %
        eval(['FR_inJC = [allFRs.',classname,'.positive.JC; allFRs.',classname,'.negative.JC];'])
        eval(['FR_inSO = [allFRs.',classname,'.positive.SO; allFRs.',classname,'.negative.SO];'])
        eval(['nonzero_inJC = [allnonzeros.',classname,'.positive.JC; allnonzeros.',classname,'.negative.JC];'])
        eval(['nonzero_inSO = [allnonzeros.',classname,'.positive.SO; allnonzeros.',classname,'.negative.SO];'])
        eval(['slope_inJC = [allslopes.',classname,'.positive.JC; -allslopes.',classname,'.negative.JC];'])
        eval(['slope_inSO = [allslopes.',classname,'.positive.SO; -allslopes.',classname,'.negative.SO];'])
        eval(['intcept_inJC = [allintcepts.',classname,'.positive.JC; allintcepts.',classname,'.negative.JC];'])
        eval(['intcept_inSO = [allintcepts.',classname,'.positive.SO; allintcepts.',classname,'.negative.SO];'])
        eval(['maxintc_inJC = [allmaxinters.',classname,'.positive.JC; allmaxinters.',classname,'.negative.JC];'])
        eval(['maxintc_inSO = [allmaxinters.',classname,'.positive.SO; allmaxinters.',classname,'.negative.SO];'])
        eval(['steepness_inJC = [allSteepness.',classname,'.positive.JC; allSteepness.',classname,'.negative.JC];'])
        eval(['steepness_inSO = [allSteepness.',classname,'.positive.SO; allSteepness.',classname,'.negative.SO];'])
        eval(['orderbias_inSO = [allOrderbias.',classname,'.positive.SO; allOrderbias.',classname,'.negative.SO];'])
        eval(['sessrange_inJC = [allsessrange.',classname,'.positive.JC; allsessrange.',classname,'.negative.JC];'])
        eval(['sessrange_inSO = [allsessrange.',classname,'.positive.SO; allsessrange.',classname,'.negative.SO];'])
        eval(['Rsq_inJC = [allRsq.',classname,'.positive.JC; allRsq.',classname,'.negative.JC];'])
        eval(['Rsq_inSO = [allRsq.',classname,'.positive.SO; allRsq.',classname,'.negative.SO];'])   
        eval(['cellnames_iclass = {allcellnames.',classname,'.positive{:}, allcellnames.',classname,'.negative{:}}'';'])
        eval(['monkeyname_list = {allmonkeys.',classname,'.positive{:}, allmonkeys.',classname,'.negative{:}}'';'])
        
        actrange_inJC = abs(maxintc_inJC - intcept_inJC);
        actrange_inSO = abs(maxintc_inSO - intcept_inSO);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        ind_G = ismember(monkeyname_list,'G');
        ind_J = ismember(monkeyname_list,'J');
        % remove outlier based on steepness; kout IQR methods
        kout = 1.5;
        % % Gerinvho
        steepness_inJC_G =  steepness_inJC(ind_G);
        quantiles_JC = quantile(steepness_inJC_G,[0.25 0.5 0.75]);
        IQR_JC = quantiles_JC(3) - quantiles_JC(1);
        steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
        steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
        %
        steepness_inSO_G =  steepness_inSO(ind_G);
        quantiles_SO = quantile(steepness_inSO_G,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];     
        steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
        %
        orderbias_inSO_G =  orderbias_inSO(ind_G);
        quantiles_SO = quantile(orderbias_inSO_G,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        orderoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];             
        %
        deltasteep = steepness_inJC - steepness_inSO;
        deltasteep_G = deltasteep(ind_G);
        quantiles_del = quantile(deltasteep_G,[0.25 0.5 0.75]);
        IQR_del = quantiles_del(3) - quantiles_del(1);
        steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
        steepoutlier_del(2) = min([steepoutlier_del(2),20]);
        %
        ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_G;
        ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_G;
        ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_G;
        ind_noneout_G = ind_goodJC & ind_goodSO; % & ind_gooddelta;       
        ind_goodorder_G = (orderbias_inSO > orderoutlier_SO(1)  & orderbias_inSO < orderoutlier_SO(2)) & ind_G;
        % % 
        % % Juan
        steepness_inJC_J =  steepness_inJC(ind_J);
        quantiles_JC = quantile(steepness_inJC_J,[0.25 0.5 0.75]);
        IQR_JC = quantiles_JC(3) - quantiles_JC(1);
        steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
        steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
        %
        steepness_inSO_J =  steepness_inSO(ind_J);
        quantiles_SO = quantile(steepness_inSO_J,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];  
        steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
        %
        orderbias_inSO_J =  orderbias_inSO(ind_J);
        quantiles_SO = quantile(orderbias_inSO_J,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        orderoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];
        %
        deltasteep = steepness_inJC - steepness_inSO;
        deltasteep_J = deltasteep(ind_J);
        quantiles_del = quantile(deltasteep_J,[0.25 0.5 0.75]);
        IQR_del = quantiles_del(3) - quantiles_del(1);
        steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
        steepoutlier_del(2) = min([steepoutlier_del(2),20]);
        %
        ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_J;
        ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_J;
        ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_J;
        ind_noneout_J = ind_goodJC & ind_goodSO; % & ind_gooddelta;  
        ind_goodorder_J = (orderbias_inSO > orderoutlier_SO(1)  & orderbias_inSO < orderoutlier_SO(2)) & ind_J;
        % %
        ind_noneout = ind_noneout_G | ind_noneout_J;
        ind_goodorder = ind_goodorder_G | ind_goodorder_J;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        %
        plot_pop_neurocorr_orderbias_JCSOinTT
    end
end


% %
% pool cell groups A and B; positive and negative
% %
if juiceABpooled & slopepooled
    cellclasses = {'OV','CV'};
%     cellclasses = {'OV'};
    ncellclasses = size(cellclasses,2);
    for icellclass = 1:ncellclasses
        classname = cellclasses{icellclass};                    
        slopesignname = 'slopemerged';
        if isequal(classname,'OV')
            FR_inJC = [allFRs.OVA.positive.JC; allFRs.OVB.positive.JC; allFRs.OVA.negative.JC; allFRs.OVB.negative.JC];
            FR_inSO = [allFRs.OVA.positive.SO; allFRs.OVB.positive.SO; allFRs.OVA.negative.SO; allFRs.OVB.negative.SO];
            nonzero_inJC = [allnonzeros.OVA.positive.JC; allnonzeros.OVB.positive.JC; allnonzeros.OVA.negative.JC; allnonzeros.OVB.negative.JC];
            nonzero_inSO = [allnonzeros.OVA.positive.SO; allnonzeros.OVB.positive.SO; allnonzeros.OVA.negative.SO; allnonzeros.OVB.negative.SO];
            slope_inJC = [allslopes.OVA.positive.JC; allslopes.OVB.positive.JC; -allslopes.OVA.negative.JC; -allslopes.OVB.negative.JC];
            slope_inSO = [allslopes.OVA.positive.SO; allslopes.OVB.positive.SO; -allslopes.OVA.negative.SO; -allslopes.OVB.negative.SO];
            intcept_inJC = [allintcepts.OVA.positive.JC; allintcepts.OVB.positive.JC; allintcepts.OVA.negative.JC; allintcepts.OVB.negative.JC];
            intcept_inSO = [allintcepts.OVA.positive.SO; allintcepts.OVB.positive.SO; allintcepts.OVA.negative.SO; allintcepts.OVB.negative.SO];
            maxintc_inJC = [allmaxinters.OVA.positive.JC; allmaxinters.OVB.positive.JC; allmaxinters.OVA.negative.JC; allmaxinters.OVB.negative.JC];
            maxintc_inSO = [allmaxinters.OVA.positive.SO; allmaxinters.OVB.positive.SO; allmaxinters.OVA.negative.SO; allmaxinters.OVB.negative.SO];
            steepness_inJC = [allSteepness.OVA.positive.JC; allSteepness.OVB.positive.JC; allSteepness.OVA.negative.JC; allSteepness.OVB.negative.JC];
            steepness_inSO = [allSteepness.OVA.positive.SO; allSteepness.OVB.positive.SO; allSteepness.OVA.negative.SO; allSteepness.OVB.negative.SO];
            orderbias_inSO = [allOrderbias.OVA.positive.SO; allOrderbias.OVB.positive.SO; allOrderbias.OVA.negative.SO; allOrderbias.OVB.negative.SO];
            sessrange_inJC = [allsessrange.OVA.positive.JC; allsessrange.OVB.positive.JC; allsessrange.OVA.negative.JC; allsessrange.OVB.negative.JC];
            sessrange_inSO = [allsessrange.OVA.positive.SO; allsessrange.OVB.positive.SO; allsessrange.OVA.negative.SO; allsessrange.OVB.negative.SO];
            Rsq_inJC = [allRsq.OVA.positive.JC; allRsq.OVB.positive.JC; allRsq.OVA.negative.JC; allRsq.OVB.negative.JC];
            Rsq_inSO = [allRsq.OVA.positive.SO; allRsq.OVB.positive.SO; allRsq.OVA.negative.SO; allRsq.OVB.negative.SO];       
            cellnames_iclass = {allcellnames.OVA.positive{:},allcellnames.OVB.positive{:},allcellnames.OVA.negative{:},allcellnames.OVB.negative{:}}';
            monkeyname_list = {allmonkeys.OVA.positive{:},allmonkeys.OVB.positive{:},allmonkeys.OVA.negative{:},allmonkeys.OVB.negative{:}}';
        elseif isequal(classname,'CV')
            FR_inJC = [allFRs.CV.positive.JC; allFRs.CV.negative.JC];
            FR_inSO = [allFRs.CV.positive.SO; allFRs.CV.negative.SO];
            nonzero_inJC = [allnonzeros.CV.positive.JC; allnonzeros.CV.negative.JC];
            nonzero_inSO = [allnonzeros.CV.positive.SO; allnonzeros.CV.negative.SO];
            slope_inJC = [allslopes.CV.positive.JC; -allslopes.CV.negative.JC];
            slope_inSO = [allslopes.CV.positive.SO; -allslopes.CV.negative.SO];
            intcept_inJC = [allintcepts.CV.positive.JC; allintcepts.CV.negative.JC];
            intcept_inSO = [allintcepts.CV.positive.SO; allintcepts.CV.negative.SO];
            maxintc_inJC = [allmaxinters.CV.positive.JC; allmaxinters.CV.negative.JC];
            maxintc_inSO = [allmaxinters.CV.positive.SO; allmaxinters.CV.negative.SO];
            steepness_inJC = [allSteepness.CV.positive.JC; allSteepness.CV.negative.JC];
            steepness_inSO = [allSteepness.CV.positive.SO; allSteepness.CV.negative.SO];
            orderbias_inSO = [allOrderbias.CV.positive.SO; allOrderbias.CV.negative.SO];
            sessrange_inJC = [allsessrange.CV.positive.JC; allsessrange.CV.negative.JC];
            sessrange_inSO = [allsessrange.CV.positive.SO; allsessrange.CV.negative.SO];
            Rsq_inJC = [allRsq.CV.positive.JC; allRsq.CV.negative.JC];
            Rsq_inSO = [allRsq.CV.positive.SO; allRsq.CV.negative.SO];       
            cellnames_iclass = {allcellnames.CV.positive{:},allcellnames.CV.negative{:}}';
            monkeyname_list = {allmonkeys.CV.positive{:},allmonkeys.CV.negative{:}}';      
            %
            slopes_QAB_inSO = [allslopes_QAB.CV.positive.SO; -allslopes_QAB.CV.negative.SO];
            slopes_chV12_inSO = [allslopes_chV12.CV.positive.SO; -allslopes_chV12.CV.negative.SO];
            intcepts_chV12_inSO = [allintcepts_chV12.CV.positive.SO; allintcepts_chV12.CV.negative.SO];
            nonzeros_QAB_inSO = [allnonzeros_QAB.CV.positive.SO; allnonzeros_QAB.CV.negative.SO];
            nonzeros_chV12_inSO = [allnonzeros_chV12.CV.positive.SO; allnonzeros_chV12.CV.negative.SO];          
        elseif isequal(classname,'CJ')
            FR_inJC = [allFRs.CJ.positive.JC; allFRs.CJ.negative.JC];
            FR_inSO = [allFRs.CJ.positive.SO; allFRs.CJ.negative.SO];
            nonzero_inJC = [allnonzeros.CJ.positive.JC; allnonzeros.CJ.negative.JC];
            nonzero_inSO = [allnonzeros.CJ.positive.SO; allnonzeros.CJ.negative.SO];
            slope_inJC = [allslopes.CJ.positive.JC; -allslopes.CJ.negative.JC];
            slope_inSO = [allslopes.CJ.positive.SO; -allslopes.CJ.negative.SO];
            slope_inSO(:,2) = -slope_inSO(:,2);
            intcept_inJC = [allintcepts.CJ.positive.JC; allintcepts.CJ.negative.JC];
            intcept_inSO = [allintcepts.CJ.positive.SO; allintcepts.CJ.negative.SO];
            maxintc_inJC = [allmaxinters.CJ.positive.JC; allmaxinters.CJ.negative.JC];
            maxintc_inSO = [allmaxinters.CJ.positive.SO; allmaxinters.CJ.negative.SO];
            steepness_inJC = [allSteepness.CJ.positive.JC; allSteepness.CJ.negative.JC];
            steepness_inSO = [allSteepness.CJ.positive.SO; allSteepness.CJ.negative.SO];
            orderbias_inSO = [allOrderbias.CJ.positive.SO; allOrderbias.CJ.negative.SO];
            sessrange_inJC = [allsessrange.CJ.positive.JC; allsessrange.CJ.negative.JC];
            sessrange_inSO = [allsessrange.CJ.positive.SO; allsessrange.CJ.negative.SO];
            Rsq_inJC = [allRsq.CJ.positive.JC; allRsq.CJ.negative.JC];
            Rsq_inSO = [allRsq.CJ.positive.SO; allRsq.CJ.negative.SO];       
            cellnames_iclass = {allcellnames.CJ.positive{:},allcellnames.CJ.negative{:}}';
            monkeyname_list = {allmonkeys.CJ.positive{:},allmonkeys.CJ.negative{:}}';      
        end
        actrange_inJC = abs(maxintc_inJC - intcept_inJC);
        actrange_inSO = abs(maxintc_inSO - intcept_inSO);
              
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % %
        ind_G = ismember(monkeyname_list,'G');
        ind_J = ismember(monkeyname_list,'J');
        % remove outlier based on steepness; kout IQR methods
        kout = 1.5;
        % % Gerinvho
        steepness_inJC_G =  steepness_inJC(ind_G);
        quantiles_JC = quantile(steepness_inJC_G,[0.25 0.5 0.75]);
        IQR_JC = quantiles_JC(3) - quantiles_JC(1);
        steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
        steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
        %
        steepness_inSO_G =  steepness_inSO(ind_G);
        quantiles_SO = quantile(steepness_inSO_G,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];     
        steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
        %
        orderbias_inSO_G =  orderbias_inSO(ind_G);
        quantiles_SO = quantile(orderbias_inSO_G,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        orderoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];             
        %
        deltasteep = steepness_inJC - steepness_inSO;
        deltasteep_G = deltasteep(ind_G);
        quantiles_del = quantile(deltasteep_G,[0.25 0.5 0.75]);
        IQR_del = quantiles_del(3) - quantiles_del(1);
        steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
        steepoutlier_del(2) = min([steepoutlier_del(2),20]);
        %
        ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_G;
        ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_G;
        ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_G;
        ind_noneout_G = ind_goodJC & ind_goodSO; % & ind_gooddelta;       
        ind_goodorder_G = (orderbias_inSO > orderoutlier_SO(1)  & orderbias_inSO < orderoutlier_SO(2)) & ind_G;
        % % 
        % % Juan
        steepness_inJC_J =  steepness_inJC(ind_J);
        quantiles_JC = quantile(steepness_inJC_J,[0.25 0.5 0.75]);
        IQR_JC = quantiles_JC(3) - quantiles_JC(1);
        steepoutlier_JC = [quantiles_JC(1)-kout*IQR_JC, quantiles_JC(3)+kout*IQR_JC];
        steepoutlier_JC(2) = min([steepoutlier_JC(2),20]);
        %
        steepness_inSO_J =  steepness_inSO(ind_J);
        quantiles_SO = quantile(steepness_inSO_J,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];  
        steepoutlier_SO(2) = min([steepoutlier_SO(2),20]);
        %
        orderbias_inSO_J =  orderbias_inSO(ind_J);
        quantiles_SO = quantile(orderbias_inSO_J,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        orderoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];
        %
        deltasteep = steepness_inJC - steepness_inSO;
        deltasteep_J = deltasteep(ind_J);
        quantiles_del = quantile(deltasteep_J,[0.25 0.5 0.75]);
        IQR_del = quantiles_del(3) - quantiles_del(1);
        steepoutlier_del = [quantiles_del(1)-kout*IQR_del, quantiles_del(3)+kout*IQR_del];
        steepoutlier_del(2) = min([steepoutlier_del(2),20]);
        %
        ind_gooddelta = (deltasteep  > steepoutlier_del(1) & deltasteep < steepoutlier_del(2)) & ind_J;
        ind_goodJC = (steepness_inJC > steepoutlier_JC(1)  & steepness_inJC < steepoutlier_JC(2)) & ind_J;
        ind_goodSO = (steepness_inSO > steepoutlier_SO(1)  & steepness_inSO < steepoutlier_SO(2)) & ind_J;
        ind_noneout_J = ind_goodJC & ind_goodSO; % & ind_gooddelta;  
        ind_goodorder_J = (orderbias_inSO > orderoutlier_SO(1)  & orderbias_inSO < orderoutlier_SO(2)) & ind_J;
        % %
        ind_noneout = ind_noneout_G | ind_noneout_J;
        ind_goodorder = ind_goodorder_G | ind_goodorder_J;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
       
        %
        plot_pop_neurocorr_orderbias_JCSOinTT
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
