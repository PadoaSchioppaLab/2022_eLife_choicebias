% pop_rho_neuronal_JCSOinTT_CVonly.m
%
% This script analyzes whether rho_neuronal encoded by chosen value cells correlate with behavior in Task 2 (SO).
% rho_neuronal in: postoffer, postjuice TW in Task 1 (JC)
%                  postoffer1 in Task 2 (SO)
%                  postoffer2 in Task 2 (SO)
%                  AB trials in Task 2 (SO)
%                  BA trials in Task 2 (SO)
%
% This script uses plot_pop_rho_neuronal_JCSOinTT.m to plot 

% author:      Jan 2021: WS

close all
clearvars

brainarea = 'OFC'; % 'DLPFC', 'VLPFC' 'OFC'
monkey_ana = 'both'; % 'Gervinho', 'Juan', 'both'

% % % 
JCSOclassifySep = 0; % if 1, load neurons that is classified seperately by JC or SO trials: 0-JCSOcellist
doSOclassifyOnly = 1;
dodoublereg = 0;

% % % 
do_bhvlogit = ''; % '': probit; '_logit'; '_neworderbias'
doshorteroffer12 ='';   % '_shorterOff12': do 250ms time windows instead of 500ms
do_CV12 = ''; % '_CV12': offer value 1 and 2 only when 1 and 2 are finally chosen

slopepooled = 1;     % if 1, pool cell postive and negative

atleast_nntrials = 2;



saveplots = 0;
if saveplots
if slopepooled
    figuresave = ['C:\Experiments\TwoTasks\Analysis\Analysis_OFC\pop_rho_neuronal_JCSOinTT_redoslope_slopepooled_',monkey_ana];
elseif ~slopepooled
    figuresave = ['C:\Experiments\TwoTasks\Analysis\Analysis_OFC\pop_rho_neuronal_JCSOinTT_redoslope_',monkey_ana];    
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

% % % % % % 
% remove outlier based on steepness across all sessions (keep consistent with behavioral measurement)
removesessions = 0; % remove based on dynamic range and saturation
dosteepout = 1; 
% Gervinho
if removesessions
    filename = ['pop_behav_ana_summary_Gervinho_removesessions',do_bhvlogit];
elseif ~removesessions
    filename = ['pop_behav_ana_summary_Gervinho',do_bhvlogit];
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
    filename = ['pop_behav_ana_summary_Juan_removesessions',do_bhvlogit];
elseif ~removesessions
    filename = ['pop_behav_ana_summary_Juan',do_bhvlogit];
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
% bhv
allSteepness  = [];
allOrderbias  = [];
allrhos       = [];
%
allmonkeys    = [];
allcellnames  = [];
% CV cells
allslopes_QAB   = [];     
allnonzeros_QAB = [];   
allintcepts_QAB = [];     

%
if JCSOclassifySep==1
    JCclasses = [JCcellist.infocell.subclass]';
    JCslopes  = [JCcellist.infocell.slopesign]';
    SOclasses = [SOcellist.infocell.subclass]';
    SOslopes  = [SOcellist.infocell.slopesign]';
    cellnames = {JCcellist.infocell.cellname}';
elseif JCSOclassifySep==0
    JCclasses = [JCSOcellist.infocell.subclass]';
    JCslopes  = [JCSOcellist.infocell.slopesign]';
    SOclasses = [JCSOcellist.infocell.subclass]';
    SOslopes  = [JCSOcellist.infocell.slopesign]';
    cellnames = {JCSOcellist.infocell.cellname}';
end

%
if     JCSOclassifySep == 0                         % classify based on both JC and SO trials - JCSOlist
    cellclassnames = {'CV'};
elseif JCSOclassifySep == 1 & doSOclassifyOnly == 0 % classify seperately based on JC and SO trials, and take the overlapped cells
    cellclassnames = {'CV_JCSOoverlap'};
elseif JCSOclassifySep == 1 & doSOclassifyOnly == 1 % classify seperately based on only SO trials
    cellclassnames = {'CV_SOonly'};    
end

cellclassnums = [3]; % CV CELLS
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
            if dodoublereg
                 filename = ['pop_rho_neuronal_JCSOinTT_redoslope',differentRho,'_',monkey_ana,'_doureg',do_bhvlogit,doshorteroffer12,do_CV12];
                % filename = ['pop_rho_neuronal_JCSOinTT_redoslope',differentRho,'_',monkey_ana,'_doureg_neworderbias'];
            else
                 filename = ['pop_rho_neuronal_JCSOinTT_redoslope',differentRho,'_',monkey_ana,do_bhvlogit,doshorteroffer12,do_CV12]; 
%                 filename = ['pop_rho_neuronal_JCSOinTT_redoslope',differentRho,'_',monkey_ana,'_neworderbias'];
            end
            load(filename);     
            %
            eval(['steepness_inJC = allSteepness.',classname,'.',slopesignname,'.JC;'])
            eval(['steepness_inSO = allSteepness.',classname,'.',slopesignname,'.SO;'])
            eval(['orderbias_inSO = allOrderbias.',classname,'.',slopesignname,'.SO;'])
            eval(['rhos_inJC = allrhos.',classname,'.',slopesignname,'.JC;'])
            eval(['rhos_inSO = allrhos.',classname,'.',slopesignname,'.SO;'])
            %
            eval(['cellnames_iclass = allcellnames.',classname,'.',slopesignname,';'])
            eval(['monkeyname_list = allmonkeys.',classname,'.',slopesignname,';'])
            %
            eval(['slopes_QAB_inJC = allslopes_QAB.',classname,'.',slopesignname,'.JC;'])
            eval(['intcepts_QAB_inJC = allintcepts_QAB.',classname,'.',slopesignname,'.JC;'])
            eval(['nonzeros_QAB_inJC = allnonzeros_QAB.',classname,'.',slopesignname,'.JC;']) 
            eval(['slopes_QAB_inSO = allslopes_QAB.',classname,'.',slopesignname,'.SO;'])
            eval(['intcepts_QAB_inSO = allintcepts_QAB.',classname,'.',slopesignname,'.SO;'])
            eval(['nonzeros_QAB_inSO = allnonzeros_QAB.',classname,'.',slopesignname,'.SO;'])   
            %
%             exnovo
        catch
            
            ind_JCiclass = JCclasses==cellclassnums(iclassname) & JCslopes==slopesigns(islopesign);
            ind_SOiclass = SOclasses==cellclassnums(iclassname) & SOslopes==slopesigns(islopesign);           
            if     JCSOclassifySep == 0                         % classify based on both JC and SO trials - JCSOlist
                ind_iclass = ind_JCiclass & ind_SOiclass;
            elseif JCSOclassifySep == 1 & doSOclassifyOnly == 0 % classify seperately based on JC and SO trials, and take the overlapped cells
                ind_iclass = ind_JCiclass & ind_SOiclass;
            elseif JCSOclassifySep == 1 & doSOclassifyOnly == 1 % classify seperately based on only SO trials
                ind_iclass = ind_SOiclass;
            end
      
            %
            ncells = sum(ind_iclass);
            cellnames_iclass = cellnames(ind_iclass);
            cellnum_iclass = find(ind_iclass==1);
            %
            monkeyname_list = [];
            % 
            steepness_inJC = []; 
            steepness_inSO = [];
            rhos_inJC = [];      
            rhos_inSO = [];           
            orderbias_inSO = [];
            %
            slopes_QAB_inJC = [];
            intcepts_QAB_inJC = [];
            nonzeros_QAB_inJC = [];
            slopes_QAB_inSO = [];
            intcepts_QAB_inSO = [];
            nonzeros_QAB_inSO = [];
            
            
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
                        rhos_inJC(icell,:) = (psyphycell.sigmoidfit.JC{3}(1)+psyphycell.sigmoidfit.SO{3}(1))/2;
                        rhos_inSO(icell,:) = (psyphycell.sigmoidfit.JC{3}(1)+psyphycell.sigmoidfit.SO{3}(1))/2;
                    elseif ~isempty(differentRho)
                        rhos_inJC(icell,:) = psyphycell.sigmoidfit.JC{3}(1);                        
                        rhos_inSO(icell,:) = psyphycell.sigmoidfit.SO{3}(1); 
                    end
                    steepness_inJC(icell,:) = psyphycell.sigmoidfit.JC{2}(2);
                    steepness_inSO(icell,:) = psyphycell.sigmoidfit.SO{2}(2);
                    orderbias_inSO(icell,:) = -2.*rhos_inSO(icell,:).*psyphycell.sigmoidfit.SO{2}(3)./psyphycell.sigmoidfit.SO{2}(2);    
                catch
                    warning off
                    [psyphycell] = sigmoidfit_TT_OrdChHyst([cellname],'logit','Both','log',1); % 'Only': SO on SO or JC on JC; 'Both': SO and JC on either SO or JC
%                     [psyphycell] = sigmoidfit_TT_OrdChHyst([cellname],'probit','Only',1); % 'Only': SO on SO or JC on JC; 'Both': SO and JC on either SO or JC         
                    relvalue_JC = exp(-(psyphycell.JC.NonChHyst.sigmoidfit.beta(1))/(psyphycell.JC.NonChHyst.sigmoidfit.beta(2)));
                    relvalue_SO = exp(-(psyphycell.SO.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SO.NonChHyst.sigmoidfit.beta(2)));                                    
                    if isempty(differentRho)
                        rhos_inJC(icell,:) = (relvalue_SO + relvalue_JC)/2;
                        rhos_inSO(icell,:) = (relvalue_SO + relvalue_JC)/2;
                    elseif ~isempty(differentRho)
                        rhos_inJC(icell,:) = relvalue_JC;
                        rhos_inSO(icell,:) = relvalue_SO;
                    end 
                    steepness_inJC(icell,:) = psyphycell.JC.NonChHyst.sigmoidfit.beta(2);
                    steepness_inSO(icell,:) = psyphycell.SO.NonChHyst.sigmoidfit.beta(2);                   
                    orderbias_inSO(icell,:) = -2.*rhos_inSO(icell,:).*psyphycell.SO.NonChHyst.sigmoidfit.beta(3)./psyphycell.SO.NonChHyst.sigmoidfit.beta(2);
                    % orderbias_inSO(icell,:) = -2.*rhos_inSO(icell,:).*psyphycell.SO.OrdChHyst.sigmoidfit.beta(3)./psyphycell.SO.OrdChHyst.sigmoidfit.beta(2);
%                     rho_inAB = exp(-(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOAB.NonChHyst.sigmoidfit.beta(2)));
%                     rho_inBA = exp(-(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(1))/(psyphycell.SOBA.NonChHyst.sigmoidfit.beta(2)));
%                     orderbias_inSO(icell,:) = rho_inBA - rho_inAB;
                end                
                %
                
                % % % 
                filename = [dirroot,cellname,'_tuning'];
                eval(['load ',filename])
                
                % redo the linear regression
                % three time windows: JC: postoffer latedelay postjuice               
                twins = {'postoffer', 'latedelay', 'postjuice'};
                ntwins = length(twins);
                for iTW = 1:ntwins
                    twin = twins{iTW};
                    eval(['neuract_JC = tuning.JC.AB.neuract.bytrial.',twin,';'])
                    act_JC = neuract_JC(:,7);
                    QVA_JC = neuract_JC(:,2);
                    QVB_JC = neuract_JC(:,3);                    
                    CJ_JC = neuract_JC(:,4).*neuract_JC(:,5);  % -1: B             
                    %
                    % average over trial type: [QVA, QVB, CJ]
                    % add other parameters into trial type table
                    % [QVA, QVB, CJ, #trials, propB, FR, std, relvalue]
                    newTrialTypes_all = [QVA_JC, QVB_JC, CJ_JC];
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
                        TrialTypes_tbl(inewTriType,8) = rhos_inJC(icell);      
                    end
                    ind_badtrialtype = TrialTypes_tbl(:,4)<atleast_nntrials;
                    TrialTypes_tbl(ind_badtrialtype,:) = [];                                                     
                    %
                    % define independent variables (use same name as cell type)                       
                    QVA = TrialTypes_tbl(:,1);
                    QVB = TrialTypes_tbl(:,2);
                    CJ = TrialTypes_tbl(:,3);   % -1 B; 1 A 
                    QVA_chA = QVA;
                    QVA_chA(CJ ==-1) = nan;
                    QVB_chB = QVB;
                    QVB_chB(CJ == 1) = nan;
                    % fitting
                    if ~dodoublereg               
                        % chosen quantity A
                        XXX = QVA_chA;
                        YYY = TrialTypes_tbl(:,6);
                        ind_good = ~isnan(XXX);
                        XXX = XXX(ind_good);
                        YYY = YYY(ind_good);
                        g = fittype('m*x+q','coeff',{'m','q'});
                        [linfit,goodness] = fit(XXX,YYY,g);
                        eval(['slopes_QAB_inJC.',twin,'(icell,1) = linfit.m;'])
                        eval(['intcepts_QAB_inJC.',twin,'(icell,1) = linfit.q;'])
                        interval = confint(linfit,0.90);
                        if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                        eval(['nonzeros_QAB_inJC.',twin,'(icell,1) = nonzero;'])
                        %
                        % chosen quantity B
                        XXX = QVB_chB;
                        YYY = TrialTypes_tbl(:,6);
                        ind_good = ~isnan(XXX);
                        XXX = XXX(ind_good);
                        YYY = YYY(ind_good);
                        g = fittype('m*x+q','coeff',{'m','q'});
                        [linfit,goodness] = fit(XXX,YYY,g);
                        eval(['slopes_QAB_inJC.',twin,'(icell,2) = linfit.m;'])
                        eval(['intcepts_QAB_inJC.',twin,'(icell,2) = linfit.q;'])
                        interval = confint(linfit,0.90);
                        if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                        eval(['nonzeros_QAB_inJC.',twin,'(icell,2) = nonzero;'])
                    else
                        XXX1 = QVA_chA;
                        XXX1(isnan(XXX1)) = 0;
                        XXX2 = QVB_chB;
                        XXX2(isnan(XXX2)) = 0;
                        YYY = TrialTypes_tbl(:,6);
%                         g = fittype('a0 + m*(k*x1 + x2)','ind',{'x1' 'x2'},'dep','y','coeff',{'a0','m','k'})
%                         [linfit,goodness] = fit([XXX1,XXX2],YYY,g);       
                        [b,bint] = regress(YYY,[ones(size(XXX1)), XXX1, XXX2],0.10);
                        eval(['slopes_QAB_inJC.',twin,'(icell,1) = b(2);'])
                        eval(['intcepts_QAB_inJC.',twin,'(icell,1) = b(1);'])
                        if prod(sign(bint(2,1)))==1, nonzero = 1; else nonzero = 0; end
                        eval(['nonzeros_QAB_inJC.',twin,'(icell,1) = nonzero;'])
                        % 
                        eval(['slopes_QAB_inJC.',twin,'(icell,2) = b(3);'])
                        eval(['intcepts_QAB_inJC.',twin,'(icell,2) = b(1);'])
                        if prod(sign(bint(3,1)))==1, nonzero = 1; else nonzero = 0; end
                        eval(['nonzeros_QAB_inJC.',twin,'(icell,2) = nonzero;'])
                    end        
                end % for itw
                
                
                % %
                % redo the linear regression
                % three time windows: SO: postoffer1 postoffer2 postjuice             
                % 
                % generate TrialTypes_tbl for each time window and define independent variables
                twins = {'postoffer1','postoffer2','postjuice'};
                QA_conds = {{'inSOoff1', 'inAB12'};{'inSOoff2', 'inBA12'};{'inSOpj','inABpj','inBApj'}};
                QB_conds = {{'inSOoff1', 'inBA12'};{'inSOoff2', 'inAB12'};{'inSOpj','inABpj','inBApj'}};
                ntwins = length(twins);
                for iTW = 1:ntwins
                    twin = twins{iTW};
                    QA_conds_iTW = QA_conds{iTW,1};
                    QB_conds_iTW = QB_conds{iTW,1};
                    nconds = length(QA_conds_iTW);
                    %
                    eval(['neuract_SO = tuning.SO.ABA.neuract.bytrial.',twin,';'])
                    act_SO = neuract_SO(:,7);
                    % do shorter time window for offer1 and offer2
                    if ~isempty(doshorteroffer12)
                        filename = [dirroot,cellname,'_data'];
                        eval(['load ',filename])
                        clear celldataerror trace trialRecord
                        ind = ismember(psyphydata(:,3),[43,44,45]);
                        psyphydata(ind,3) = 50;
                        %
                        if isequal(twin,'postoffer1')
                             timewindow_new = {'offer1_shorter',[30 30],[100 450]};
                             act_SO = actmake_TT(celldata, psyphydata, neuract_SO(:,1), timewindow_new);
                        elseif isequal(twin,'postoffer2')
                            timewindow_new = {'offer2_shorter',[32 32],[100 450]};
                             act_SO = actmake_TT(celldata, psyphydata, neuract_SO(:,1), timewindow_new);
                        end
                    end                 
                    %                   
                    QVA_SO = neuract_SO(:,2);
                    QVB_SO = neuract_SO(:,3);                    
                    CJ_SO = neuract_SO(:,5);  % -1: B
                    ord_SO = neuract_SO(:,6); % -1: BA
                    CO_SO = CJ_SO.*ord_SO;    % -1: chosen 2nd                               
                    %
                    % average over trial type: [QVA, QVB, ord, CJ]
                    % add other parameters into trial type table
                    % [QVA, QVB, ord, CJ, #trials, propB, FR, std, relvalue]
                    newTrialTypes_all = [QVA_SO, QVB_SO, ord_SO, CJ_SO];
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
                        TrialTypes_tbl(inewTriType,9) = rhos_inSO(icell);      
                    end
                    ind_badtrialtype = TrialTypes_tbl(:,5)<atleast_nntrials;
                    TrialTypes_tbl(ind_badtrialtype,:) = [];
                    %
                    % define independent variables (use same name as cell type)  
                    ABBA = TrialTypes_tbl(:,3); % -1 BA; 1 AB  
                    QVA = TrialTypes_tbl(:,1);
                    QVAinBA = QVA;
                    QVAinBA(ABBA== 1) = nan;
                    QVAinAB = QVA;
                    QVAinAB(ABBA==-1) = nan; 
                    QVB = TrialTypes_tbl(:,2);
                    QVBinBA = QVB;
                    QVBinBA(ABBA== 1) = nan;
                    QVBinAB = QVB;
                    QVBinAB(ABBA==-1) = nan;                        
                    chJ = TrialTypes_tbl(:,4); % -1 B; 1 A  
                    QVA_chA = QVA;
                    QVA_chA(chJ ==-1) = nan;
                    QVA_chA_inAB = QVA_chA;
                    QVA_chA_inAB(ABBA ==-1) = nan;
                    QVA_chA_inBA = QVA_chA;
                    QVA_chA_inBA(ABBA == 1) = nan;                    
                    QVB_chB = QVB;
                    QVB_chB(chJ == 1) = nan;
                    QVB_chB_inAB = QVB_chB;
                    QVB_chB_inAB(ABBA ==-1) = nan;
                    QVB_chB_inBA = QVB_chB;
                    QVB_chB_inBA(ABBA == 1) = nan;
                    %
                    if iTW == 1
                        if isempty(do_CV12)
                            QVA_inSOoff1 = QVAinAB;
                            QVA_inAB12   = QVAinAB;                       
                            QVB_inSOoff1 = QVBinBA; 
                            QVB_inBA12   = QVBinBA; 
                        elseif ~isempty(do_CV12)
                            QVA_inSOoff1 = QVA_chA_inAB;
                            QVA_inAB12   = QVA_chA_inAB;                       
                            QVB_inSOoff1 = QVB_chB_inBA; 
                            QVB_inBA12   = QVB_chB_inBA; 
                        end
                    elseif iTW == 2 
                        if isempty(do_CV12)
                            QVA_inSOoff2 = QVAinBA;
                            QVA_inBA12   = QVAinBA;
                            QVB_inSOoff2 = QVBinAB;  
                            QVB_inAB12   = QVBinAB; 
                        elseif ~isempty(do_CV12)
                            QVA_inSOoff2 = QVA_chA_inBA;
                            QVA_inBA12   = QVA_chA_inBA;
                            QVB_inSOoff2 = QVB_chB_inAB;  
                            QVB_inAB12   = QVB_chB_inAB; 
                        end                      
                    elseif iTW == 3 
                        QVA_inSOpj = QVA_chA;
                        QVB_inSOpj = QVB_chB;
                        QVA_inABpj = QVA_chA_inAB;
                        QVB_inABpj = QVB_chB_inAB;
                        QVA_inBApj = QVA_chA_inBA;
                        QVB_inBApj = QVB_chB_inBA;
                    end
                    % fitting
                    for icond = 1:nconds
                        if ~dodoublereg             
                            % chosen quantity A
                            QA_icond = QA_conds_iTW{icond};
                            eval(['XXX = QVA_',QA_icond,';'])                   
                            YYY = TrialTypes_tbl(:,7);
                            ind_good = ~isnan(XXX);
                            XXX = XXX(ind_good);
                            YYY = YYY(ind_good);
                            % XXX(isnan(XXX)) = 0;
                            g = fittype('m*x+q','coeff',{'m','q'});
                            [linfit,goodness] = fit(XXX,YYY,g);
                            eval(['slopes_QAB_inSO.',QA_icond,'(icell,1) = linfit.m;'])
                            eval(['intcepts_QAB_inSO.',QA_icond,'(icell,1) = linfit.q;'])
                            interval = confint(linfit,0.90);
                            if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                            eval(['nonzeros_QAB_inSO.',QA_icond,'(icell,1) = nonzero;'])
                            % chosen quantity B
                            QB_icond = QB_conds_iTW{icond};
                            eval(['XXX = QVB_',QB_icond,';'])                   
                            YYY = TrialTypes_tbl(:,7);
                            ind_good = ~isnan(XXX);
                            XXX = XXX(ind_good);
                            YYY = YYY(ind_good);
                            % XXX(isnan(XXX)) = 0;
                            g = fittype('m*x+q','coeff',{'m','q'});
                            [linfit,goodness] = fit(XXX,YYY,g);
                            eval(['slopes_QAB_inSO.',QB_icond,'(icell,2) = linfit.m;'])
                            eval(['intcepts_QAB_inSO.',QB_icond,'(icell,2) = linfit.q;'])
                            interval = confint(linfit,0.90);
                            if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                            eval(['nonzeros_QAB_inSO.',QB_icond,'(icell,2) = nonzero;'])
                        elseif dodoublereg
                            %
                            QA_icond = QA_conds_iTW{icond};
                            eval(['XXX1 = QVA_',QA_icond,';'])       
                            XXX1(isnan(XXX1)) = 0;
                            QB_icond = QB_conds_iTW{icond};
                            eval(['XXX2 = QVB_',QB_icond,';'])  
                            XXX2(isnan(XXX2)) = 0;
                            YYY = TrialTypes_tbl(:,6);
    %                         g = fittype('a0 + m*(k*x1 + x2)','ind',{'x1' 'x2'},'dep','y','coeff',{'a0','m','k'})
    %                         [linfit,goodness] = fit([XXX1,XXX2],YYY,g);       
                            [b,bint] = regress(YYY,[ones(size(XXX1)), XXX1, XXX2],0.10);
                            eval(['slopes_QAB_inSO.',QB_icond,'(icell,1) = b(2);'])
                            eval(['intcepts_QAB_inSO.',QB_icond,'(icell,1) = b(1);'])
                            if prod(sign(bint(2,1)))==1, nonzero = 1; else nonzero = 0; end
                            eval(['nonzeros_QAB_inSO.',QB_icond,'(icell,1) = nonzero;'])
                            % 
                            eval(['slopes_QAB_inSO.',QB_icond,'(icell,2) = b(3);'])
                            eval(['intcepts_QAB_inSO.',QB_icond,'(icell,2) = b(1);'])
                            if prod(sign(bint(3,1)))==1, nonzero = 1; else nonzero = 0; end
                            eval(['nonzeros_QAB_inSO.',QB_icond,'(icell,2) = nonzero;'])
                            
                        end
                    end % for icond
                end % for iTW                              
            end %for icell
            %
            
            % PUT DATA TOGETHER
            eval(['allcellnames.',classname,'.',slopesignname,' = cellnames_iclass;'])
            eval(['allmonkeys.',classname,'.',slopesignname,' = monkeyname_list;'])
            %
            eval(['allSteepness.',classname,'.',slopesignname,'.JC = steepness_inJC;'])
            eval(['allSteepness.',classname,'.',slopesignname,'.SO = steepness_inSO;'])
            eval(['allOrderbias.',classname,'.',slopesignname,'.SO = orderbias_inSO;'])
            eval(['allrhos.',classname,'.',slopesignname,'.JC = rhos_inJC;'])
            eval(['allrhos.',classname,'.',slopesignname,'.SO = rhos_inSO;'])
            %
            eval(['allslopes_QAB.',classname,'.',slopesignname,'.JC = slopes_QAB_inJC;']) 
            eval(['allintcepts_QAB.',classname,'.',slopesignname,'.JC = intcepts_QAB_inJC;'])   
            eval(['allnonzeros_QAB.',classname,'.',slopesignname,'.JC = nonzeros_QAB_inJC;'])             
            eval(['allslopes_QAB.',classname,'.',slopesignname,'.SO = slopes_QAB_inSO;']) 
            eval(['allintcepts_QAB.',classname,'.',slopesignname,'.SO = intcepts_QAB_inSO;'])   
            eval(['allnonzeros_QAB.',classname,'.',slopesignname,'.SO = nonzeros_QAB_inSO;'])     
            
            % %
            if dodoublereg
                filename = ['pop_rho_neuronal_JCSOinTT_redoslope',differentRho,'_',monkey_ana,'_doureg',do_bhvlogit,doshorteroffer12,do_CV12];
                % filename = ['pop_rho_neuronal_JCSOinTT_redoslope',differentRho,'_',monkey_ana,'_doureg_neworderbias'];
            else
                filename = ['pop_rho_neuronal_JCSOinTT_redoslope',differentRho,'_',monkey_ana,do_bhvlogit,doshorteroffer12,do_CV12]; 
                % filename = ['pop_rho_neuronal_JCSOinTT_redoslope',differentRho,'_',monkey_ana,'_neworderbias'];
            end
            eval(['save ',filename ' allSteepness allOrderbias allrhos allcellnames allmonkeys  '...
                                   ' allslopes_QAB allnonzeros_QAB allintcepts_QAB '])
                   
        end % try catch end
        
        % % % 
        % plot each cell class
        % % %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        if ~slopepooled
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
        steepoutlier_JC(2) = min([steepoutlier_JC(2)]);
        %
        steepness_inSO_G =  steepness_inSO(ind_G);
        quantiles_SO = quantile(steepness_inSO_G,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];     
        steepoutlier_SO(2) = min([steepoutlier_SO(2)]);
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
        steepoutlier_del(2) = min([steepoutlier_del(2)]);
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
        steepoutlier_JC(2) = min([steepoutlier_JC(2)]);
        %
        steepness_inSO_J =  steepness_inSO(ind_J);
        quantiles_SO = quantile(steepness_inSO_J,[0.25 0.5 0.75]);
        IQR_SO = quantiles_SO(3) - quantiles_SO(1);
        steepoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];  
        steepoutlier_SO(2) = min([steepoutlier_SO(2)]);
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
        steepoutlier_del(2) = min([steepoutlier_del(2)]);
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
         
        plot_pop_rho_neuronal_JCSOinTT
        
        end
    end %for islopesign
end %for iclassname

                 
% %
% pool cell groups positive and negative
% %
if slopepooled
    %
    if     JCSOclassifySep == 0                         % classify based on both JC and SO trials - JCSOlist
        cellclasses = {'CV'};
    elseif JCSOclassifySep == 1 & doSOclassifyOnly == 0 % classify seperately based on JC and SO trials, and take the overlapped cells
        cellclasses = {'CV_JCSOoverlap'};
    elseif JCSOclassifySep == 1 & doSOclassifyOnly == 1 % classify seperately based on only SO trials
        cellclasses = {'CV_SOonly'};    
    end
    
    ncellclasses = size(cellclasses,2);
    for icellclass = 1:ncellclasses
        classname = cellclasses{icellclass};                    
        slopesignname = 'slopemerged';
        %
        eval(['steepness_inJC = [allSteepness.',classname,'.positive.JC; allSteepness.',classname,'.negative.JC];'])
        eval(['steepness_inSO = [allSteepness.',classname,'.positive.SO; allSteepness.',classname,'.negative.SO];'])
        eval(['orderbias_inSO = [allOrderbias.',classname,'.positive.SO; allOrderbias.',classname,'.negative.SO];'])
        eval(['rhos_inJC = [allrhos.',classname,'.positive.JC; allrhos.',classname,'.negative.JC];'])
        eval(['rhos_inSO = [allrhos.',classname,'.positive.SO; allrhos.',classname,'.negative.SO];'])
        %
        eval(['cellnames_iclass = {allcellnames.',classname,'.positive{:},allcellnames.',classname,'.negative{:}}'';'])
        eval(['monkeyname_list = {allmonkeys.',classname,'.positive{:},allmonkeys.',classname,'.negative{:}}'';'])
        %
        eval(['slopes_QAB_inJC.postoffer = [allslopes_QAB.',classname,'.positive.JC.postoffer; -allslopes_QAB.',classname,'.negative.JC.postoffer];'])
        eval(['slopes_QAB_inJC.latedelay = [allslopes_QAB.',classname,'.positive.JC.latedelay; -allslopes_QAB.',classname,'.negative.JC.latedelay];'])
        eval(['slopes_QAB_inJC.postjuice = [allslopes_QAB.',classname,'.positive.JC.postjuice; -allslopes_QAB.',classname,'.negative.JC.postjuice];'])
        eval(['nonzeros_QAB_inJC.postoffer = [allnonzeros_QAB.',classname,'.positive.JC.postoffer; allnonzeros_QAB.',classname,'.negative.JC.postoffer];'])  
        eval(['nonzeros_QAB_inJC.latedelay = [allnonzeros_QAB.',classname,'.positive.JC.latedelay; allnonzeros_QAB.',classname,'.negative.JC.latedelay];'])
        eval(['nonzeros_QAB_inJC.postjuice = [allnonzeros_QAB.',classname,'.positive.JC.postjuice; allnonzeros_QAB.',classname,'.negative.JC.postjuice];'])
        eval(['intcepts_QAB_inJC.postoffer = [allintcepts_QAB.',classname,'.positive.JC.postoffer; allintcepts_QAB.',classname,'.negative.JC.postoffer];'])
        eval(['intcepts_QAB_inJC.latedelay = [allintcepts_QAB.',classname,'.positive.JC.latedelay; allintcepts_QAB.',classname,'.negative.JC.latedelay];'])
        eval(['intcepts_QAB_inJC.postjuice = [allintcepts_QAB.',classname,'.positive.JC.postjuice; allintcepts_QAB.',classname,'.negative.JC.postjuice];'])
        %
        eval(['slopes_QAB_inSO.inSOoff1 = [allslopes_QAB.',classname,'.positive.SO.inSOoff1; -allslopes_QAB.',classname,'.negative.SO.inSOoff1];'])          
        eval(['slopes_QAB_inSO.inSOoff2 = [allslopes_QAB.',classname,'.positive.SO.inSOoff2; -allslopes_QAB.',classname,'.negative.SO.inSOoff2];'])
        eval(['slopes_QAB_inSO.inAB12   = [allslopes_QAB.',classname,'.positive.SO.inAB12;   -allslopes_QAB.',classname,'.negative.SO.inAB12];'])
        eval(['slopes_QAB_inSO.inBA12   = [allslopes_QAB.',classname,'.positive.SO.inBA12;   -allslopes_QAB.',classname,'.negative.SO.inBA12];'])
        eval(['slopes_QAB_inSO.inABpj   = [allslopes_QAB.',classname,'.positive.SO.inABpj;   -allslopes_QAB.',classname,'.negative.SO.inABpj];'])
        eval(['slopes_QAB_inSO.inBApj   = [allslopes_QAB.',classname,'.positive.SO.inBApj;   -allslopes_QAB.',classname,'.negative.SO.inBApj];'])
        eval(['nonzeros_QAB_inSO.inSOoff1 = [allnonzeros_QAB.',classname,'.positive.SO.inSOoff1; allnonzeros_QAB.',classname,'.negative.SO.inSOoff1];'])          
        eval(['nonzeros_QAB_inSO.inSOoff2 = [allnonzeros_QAB.',classname,'.positive.SO.inSOoff2; allnonzeros_QAB.',classname,'.negative.SO.inSOoff2];'])
        eval(['nonzeros_QAB_inSO.inAB12   = [allnonzeros_QAB.',classname,'.positive.SO.inAB12;   allnonzeros_QAB.',classname,'.negative.SO.inAB12];'])
        eval(['nonzeros_QAB_inSO.inBA12   = [allnonzeros_QAB.',classname,'.positive.SO.inBA12;   allnonzeros_QAB.',classname,'.negative.SO.inBA12];'])
        eval(['nonzeros_QAB_inSO.inABpj   = [allnonzeros_QAB.',classname,'.positive.SO.inABpj;   allnonzeros_QAB.',classname,'.negative.SO.inABpj];'])
        eval(['nonzeros_QAB_inSO.inBApj   = [allnonzeros_QAB.',classname,'.positive.SO.inBApj;   allnonzeros_QAB.',classname,'.negative.SO.inBApj];'])
        eval(['intcepts_QAB_inSO.inSOoff1 = [allintcepts_QAB.',classname,'.positive.SO.inSOoff1; allintcepts_QAB.',classname,'.negative.SO.inSOoff1]; '])         
        eval(['intcepts_QAB_inSO.inSOoff2 = [allintcepts_QAB.',classname,'.positive.SO.inSOoff2; allintcepts_QAB.',classname,'.negative.SO.inSOoff2];'])
        eval(['intcepts_QAB_inSO.inAB12   = [allintcepts_QAB.',classname,'.positive.SO.inAB12;   allintcepts_QAB.',classname,'.negative.SO.inAB12];'])
        eval(['intcepts_QAB_inSO.inBA12   = [allintcepts_QAB.',classname,'.positive.SO.inBA12;   allintcepts_QAB.',classname,'.negative.SO.inBA12];'])
        eval(['intcepts_QAB_inSO.inABpj   = [allintcepts_QAB.',classname,'.positive.SO.inABpj;   allintcepts_QAB.',classname,'.negative.SO.inABpj];'])
        eval(['intcepts_QAB_inSO.inBApj   = [allintcepts_QAB.',classname,'.positive.SO.inBApj;   allintcepts_QAB.',classname,'.negative.SO.inBApj];'])                    
        
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
%         orderbias_inSO_G =  orderbias_inSO(ind_G);
%         quantiles_SO = quantile(orderbias_inSO_G,[0.25 0.5 0.75]);
%         IQR_SO = quantiles_SO(3) - quantiles_SO(1);
%         orderoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];             
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
%         ind_goodorder_G = (orderbias_inSO > orderoutlier_SO(1)  & orderbias_inSO < orderoutlier_SO(2)) & ind_G;
        ind_goodJC = (steepness_inJC > steepoutlier_JC_G(1)  & steepness_inJC < steepoutlier_JC_G(2)) & ind_G;
        ind_goodSO = (steepness_inSO > steepoutlier_SO_G(1)  & steepness_inSO < steepoutlier_SO_G(2)) & ind_G;
        ind_noneout_G = ind_goodJC & ind_goodSO; % & ind_gooddelta;
        % % 
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
%         orderbias_inSO_J =  orderbias_inSO(ind_J);
%         quantiles_SO = quantile(orderbias_inSO_J,[0.25 0.5 0.75]);
%         IQR_SO = quantiles_SO(3) - quantiles_SO(1);
%         orderoutlier_SO = [quantiles_SO(1)-kout*IQR_SO, quantiles_SO(3)+kout*IQR_SO];
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
%         ind_noneout_J = ind_goodJC & ind_goodSO; % & ind_gooddelta;  
%         ind_goodorder_J = (orderbias_inSO > orderoutlier_SO(1)  & orderbias_inSO < orderoutlier_SO(2)) & ind_J;
        ind_goodJC = (steepness_inJC > steepoutlier_JC_J(1)  & steepness_inJC < steepoutlier_JC_J(2)) & ind_J;
        ind_goodSO = (steepness_inSO > steepoutlier_SO_J(1)  & steepness_inSO < steepoutlier_SO_J(2)) & ind_J;
        ind_noneout_J = ind_goodJC & ind_goodSO; % & ind_gooddelta;
        % %
        ind_noneout = ind_noneout_G | ind_noneout_J;
        % ind_goodorder = ind_goodorder_G | ind_goodorder_J;
        if ~dosteepout
            ind_noneout = logical(ones(size(ind_noneout)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
       
        %
        plot_pop_rho_neuronal_JCSOinTT
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
