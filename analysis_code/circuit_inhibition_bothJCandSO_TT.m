% circuit_inhibition_new_bothJCandSO_TT.m
%
% This script tests circuit inhibition in a new way, the main result is
% expected to be the same as circuit_inhibition_TT.m. 
% What are new in this script:
% 1. test the anti-corr between FR in interoffer delay and V(O)
% 2. test the circuit inhibition in JC (maybe wrong)
% 3. correlate the anti-corr and steepness 
% 4. correlate the anti-corr and order bias
% 5. correlate JC and SO
% 6. ...

% author:   Feb 2020 - WS
% revision: Feb 2021 - WS


close all
clearvars

brainarea = 'OFC'; % 'DLPFC', 'VLPFC' 'OFC'
monkey_ana = 'both'; % 'Gervinho', 'Juan', 'both'
monkey_plot = 'both';

mintrialnum = 100;

dobhvlogit = ''; % '': probit; '_logit'; '_neworderbias'

differentRho = '_differRho';    % '_differRho' or ''; do analysis based on the same rho in JC and SO or not 
twoTWinJCandSO = '';  % '' or ''; do fewer TW: 2 for JC and 2 for SO
doJCseq = '';          % '_JCseq' or ''; do sequential JC parameters or not
nntrials = '_2nntrials'; % '_2nntrials' or '_3nntrials'
anovasetup = '_01pvallessTW_95nonzero'; % '_05pval_95nonzero'; '_01pvallessTW_95nonzero';
filename = ['TTcellist',doJCseq,'_',brainarea,'_both',differentRho,nntrials,anovasetup];
load(filename);

% cell analysis types - identify cells based on different conditions
JCSOclassifySep = 0; % if 1, load neurons that is classified seperately by JC or SO trials: 0-JCSOcellist
                        % only use 'SO' in this script

JCSOs = {'SO'};
nJCSO = size(JCSOs,2);       

alignments = {'off'}; % 'off', 'cho', 'cue'
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

%
for iclass = 1:nclassnames % nclassnames == 1; CJ only
    classnamenum = classnamenums(iclass);
    ind_JCiclass = JCclasses==classnamenum;
    ind_SOiclass = SOclasses==classnamenum;
    % ind_iclass = ind_JCiclass & ind_SOiclass; % |: pass the classification in either JC or SO; &: pass the classification in both JC and SO
    ind_iclass = ind_SOiclass;
    %
    ncells = sum(ind_iclass);
    cellnames_iclass = cellnames(ind_iclass);
    cellnum_iclass   = find(ind_iclass==1);
    JCslopes_iclass  = JCslopes(ind_iclass); % 1: CJA, -1 CJB
    SOslopes_iclass  = SOslopes(ind_iclass); % 1: CJA, -1 CJB
    ncells_ana = 0;
    %
    orderbias_all_SO = [];
    Steepness_all_sumJCSO = [];
    averCritFR_vs_VO_RR_sumJCSO = [];
    averCritFR_vs_VO_PP_sumJCSO = []; 
    averCritFR_vs_VO_a1_sumJCSO = [];
    
    for iJCSO = 1:nJCSO
        JCorSO = JCSOs{iJCSO};
        [~, ~,  timewindows_off, timewindows_cho, ~, timewindows_cue, ~] = get_timewindows_TT(JCorSO);
    
        for ialign = 1:naligns
            alignment = alignments{ialign};         
            eval(['ntwins = size(timewindows_',alignment,',2);'])
            ncells_ana = 0;
            %
            timepoints = [];
            %
            FRvsVO_RR_movTW = [];
            FRvsVO_PP_movTW = [];
            %
            FRvsVO_a0_movTW = [];
            FRvsVO_a1_movTW = [];
            %
            FRVO_movTW = [];
            FRVE_movTW = [];
            %
            FRVO_control_movTW = [];
            FRVE_control_movTW = [];
            %
            a1vsOrderbias_RR_movTW = [];
            PPvsSteep_PP_movTW = [];
            %
            Steepness_all = [];    
            orderbias_all = [];
            rho_all = [];
            valuerange_all = [];
            %
            averCritFR_vs_VO_RR = []; 
            averCritFR_vs_VO_PP = [];
            %
            averCritFR_vs_VO_a0 = [];
            averCritFR_vs_VO_a1 = [];
            %
            corr_VOandVE = [];
            
            try
%                 exnovo
                if ~JCSOclassifySep
                    filename = ['pop_circuit_inhibition_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,dobhvlogit];
                    % filename = ['pop_circuit_inhibition_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,'_neworderbias'];
                elseif JCSOclassifySep
                    filename = ['pop_circuit_inhibition_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,'_SOclassonly',dobhvlogit];
%                     filename = ['pop_circuit_inhibition_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,'_SOclassonly_neworderbias'];      
                end
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
                    try % load behavior
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
                            orderbias_SO = rho_inBA - rho_inAB;
                        else
                            orderbias_SO = -2.*rho_SO.*psyphycell.SO.NonChHyst.sigmoidfit.beta(3)./psyphycell.SO.NonChHyst.sigmoidfit.beta(2);
                            % orderbias_SO = -2.*rho_SO.*psyphycell.SO.OrdChHyst.sigmoidfit.beta(3)./psyphycell.SO.OrdChHyst.sigmoidfit.beta(2);
                        end
                    end                             

                    % JC trials
                    if isequal(JCorSO,'JC')
                        pairname = 'AB';
                        eval(['superneuract=tuning.',JCorSO,'.',pairname,'.superneuract.bytrial;']) 
                        Steepness_all(:,ncells_ana) = steepness_JC;
                        rho_all(:,ncells_ana) = rho_JC;                        
                    % SO trials
                    elseif isequal(JCorSO,'SO')
                        pairname = 'ABA';
                        eval(['superneuract=tuning.',JCorSO,'.',pairname,'.superneuract.bytrial;'])
                        Steepness_all(:,ncells_ana) = steepness_SO;
                        rho_all(:,ncells_ana) = rho_SO;
                    end % if JC or SO
                    orderbias_all(:,ncells_ana) = orderbias_SO;
                                        
                    %
                    allneuroFR_OE_icell = [];
                    %
                    for itwin = 1:ntwins
                        eval(['timepoints(itwin,1) = mean(timewindows_',alignment,'{itwin}{3});'])
                        %
                        eval(['neuract = superneuract.twds_',alignment,num2str(itwin),';'])     
                        % remove force choice
                        ind_forced = neuract(:,2)==0 | neuract(:,3)==0;
                        neuract_all = neuract;
                        neuract(ind_forced,:) = [];

                        if isequal(JCorSO,'SO')
                            orderABBA = neuract(:,6);
                           if SOslope_icell == 1 % CJA
                                ind_BA = orderABBA == -1;
                                neuract_OE = neuract(ind_BA,:);
                                ind_AB = orderABBA == 1;
                                neuract_EO = neuract(ind_AB,:);
                            elseif SOslope_icell == -1 %CJB
                                ind_AB = orderABBA == 1;
                                neuract_OE = neuract(ind_AB,:);
                                ind_BA = orderABBA == -1;
                                neuract_EO = neuract(ind_BA,:);
                            end   
                        end
                        %
                        % JC trials
                        if isequal(JCorSO,'JC')
                            ChosenJuice = neuract(:,4).*neuract(:,5); % 1:chosen A; -1:chosen B                                 
                            neuract = neuract(ChosenJuice == JCslope_icell,:);
                            OfferValueA = abs(neuract(:,2)).*rho_JC;
                            OfferValueB = abs(neuract(:,3));
                            neurFR = neuract(:,7);
                            %
                            % calculate correlation between V(O) and neural FR
                            if JCslope_icell == 1 % CJA
                                [rr,pp] = corr(OfferValueB,neurFR);
                                aa = polyfit(OfferValueB,neurFR,1);
                                a0 = aa(2);
                                a1 = aa(1);                            
                            elseif JCslope_icell == -1 %CJB
                                [rr,pp] = corr(OfferValueA,neurFR);
                                aa = polyfit(OfferValueA,neurFR,1);
                                a0 = aa(2);
                                a1 = aa(1); 
                            end
                            FRvsVO_RR_movTW(itwin,ncells_ana) = rr;
                            FRvsVO_PP_movTW(itwin,ncells_ana) = pp;
                            FRvsVO_a0_movTW(itwin,ncells_ana) = a0;
                            FRvsVO_a1_movTW(itwin,ncells_ana) = a1;
                            %
                            % separate trials into three quantile V(O) high mid low
                            if JCslope_icell == 1 % CJA
                                thres = quantile(OfferValueB,2);
                                ind_low  = OfferValueB< thres(1);
                                ind_mid  = OfferValueB< thres(2) & OfferValueB>=thres(1);
                                ind_high = OfferValueB>=thres(2);
                                FRVO_movTW.low(itwin,ncells_ana)  = nanmean(neurFR(ind_low));
                                FRVO_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR(ind_mid));
                                FRVO_movTW.high(itwin,ncells_ana) = nanmean(neurFR(ind_high));   
                            elseif JCslope_icell == -1 %CJB
                                thres = quantile(OfferValueA,2);
                                ind_low  = OfferValueA< thres(1);
                                ind_mid  = OfferValueA< thres(2) & OfferValueA>=thres(1);
                                ind_high = OfferValueA>=thres(2);
                                FRVO_movTW.low(itwin,ncells_ana)  = nanmean(neurFR(ind_low));
                                FRVO_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR(ind_mid));
                                FRVO_movTW.high(itwin,ncells_ana) = nanmean(neurFR(ind_high)); 
                            end    
                        %       
                        % SO trials   
                        elseif isequal(JCorSO,'SO')
                            OfferValueA_OE = abs(neuract_OE(:,2)).*rho_SO;
                            OfferValueB_OE = abs(neuract_OE(:,3));
                            neurFR_OE = neuract_OE(:,7);
                            OfferValueA_EO = abs(neuract_EO(:,2)).*rho_SO;
                            OfferValueB_EO = abs(neuract_EO(:,3));
                            neurFR_EO = neuract_EO(:,7);
                            % calculate correlation between V(O) and neural FR
                            if SOslope_icell == 1 % CJA
                                [rr,pp] = corr(OfferValueB_OE,neurFR_OE);
                                aa = polyfit(OfferValueB_OE,neurFR_OE,1);
                                a0 = aa(2);
                                a1 = aa(1); 
                            elseif SOslope_icell == -1 %CJB
                                [rr,pp] = corr(OfferValueA_OE,neurFR_OE);
                                aa = polyfit(OfferValueA_OE,neurFR_OE,1);
                                a0 = aa(2);
                                a1 = aa(1); 
                            end
                            FRvsVO_RR_movTW(itwin,ncells_ana) = rr;
                            FRvsVO_PP_movTW(itwin,ncells_ana) = pp;
                            FRvsVO_a0_movTW(itwin,ncells_ana) = a0;
                            FRvsVO_a1_movTW(itwin,ncells_ana) = a1;
                            %
                            % separate OE trials into three quantile V(O) high mid low
                            if SOslope_icell == 1 % CJA
                                thres = quantile(OfferValueB_OE,2);
                                ind_low  = OfferValueB_OE< thres(1);
                                ind_mid  = OfferValueB_OE< thres(2) & OfferValueB_OE>=thres(1);
                                ind_high = OfferValueB_OE>=thres(2);
                                FRVO_movTW.low(itwin,ncells_ana)  = nanmean(neurFR_OE(ind_low));
                                FRVO_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR_OE(ind_mid));
                                FRVO_movTW.high(itwin,ncells_ana) = nanmean(neurFR_OE(ind_high));   
                            elseif SOslope_icell == -1 %CJB
                                thres = quantile(OfferValueA_OE,2);
                                ind_low  = OfferValueA_OE< thres(1);
                                ind_mid  = OfferValueA_OE< thres(2) & OfferValueA_OE>=thres(1);
                                ind_high = OfferValueA_OE>=thres(2);
                                FRVO_movTW.low(itwin,ncells_ana)  = nanmean(neurFR_OE(ind_low));
                                FRVO_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR_OE(ind_mid));
                                FRVO_movTW.high(itwin,ncells_ana) = nanmean(neurFR_OE(ind_high)); 
                            end   
                            % A control separation 
                            % separate OE trials into three quantile V(E) high mid low
                            if SOslope_icell == 1 % CJA
                                thres = quantile(OfferValueA_OE,2);
                                ind_low  = OfferValueA_OE< thres(1);
                                ind_mid  = OfferValueA_OE< thres(2) & OfferValueA_OE>=thres(1);
                                ind_high = OfferValueA_OE>=thres(2);
                                FRVE_control_movTW.low(itwin,ncells_ana)  = nanmean(neurFR_OE(ind_low));
                                FRVE_control_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR_OE(ind_mid));
                                FRVE_control_movTW.high(itwin,ncells_ana) = nanmean(neurFR_OE(ind_high));   
                            elseif SOslope_icell == -1 %CJB
                                thres = quantile(OfferValueB_OE,2);
                                ind_low  = OfferValueB_OE< thres(1);
                                ind_mid  = OfferValueB_OE< thres(2) & OfferValueB_OE>=thres(1);
                                ind_high = OfferValueB_OE>=thres(2);
                                FRVE_control_movTW.low(itwin,ncells_ana)  = nanmean(neurFR_OE(ind_low));
                                FRVE_control_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR_OE(ind_mid));
                                FRVE_control_movTW.high(itwin,ncells_ana) = nanmean(neurFR_OE(ind_high)); 
                            end   
                            %
                            % separate EO trials into three quantile V(E) high mid low
                            if SOslope_icell == 1 % CJA
                                thres = quantile(OfferValueA_EO,2);
                                ind_low  = OfferValueA_EO< thres(1);
                                ind_mid  = OfferValueA_EO< thres(2) & OfferValueA_EO>=thres(1);
                                ind_high = OfferValueA_EO>=thres(2);
                                FRVE_movTW.low(itwin,ncells_ana)  = nanmean(neurFR_EO(ind_low));
                                FRVE_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR_EO(ind_mid));
                                FRVE_movTW.high(itwin,ncells_ana) = nanmean(neurFR_EO(ind_high));   
                            elseif SOslope_icell == -1 %CJB
                                thres = quantile(OfferValueB_EO,2);
                                ind_low  = OfferValueB_EO< thres(1);
                                ind_mid  = OfferValueB_EO< thres(2) & OfferValueB_EO>=thres(1);
                                ind_high = OfferValueB_EO>=thres(2);
                                FRVE_movTW.low(itwin,ncells_ana)  = nanmean(neurFR_EO(ind_low));
                                FRVE_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR_EO(ind_mid));
                                FRVE_movTW.high(itwin,ncells_ana) = nanmean(neurFR_EO(ind_high)); 
                            end      
                            % A control separation 
                            % separate EO trials into three quantile V(O) high mid low
                            if SOslope_icell == 1 % CJA
                                thres = quantile(OfferValueB_EO,2);
                                ind_low  = OfferValueB_EO< thres(1);
                                ind_mid  = OfferValueB_EO< thres(2) & OfferValueB_EO>=thres(1);
                                ind_high = OfferValueB_EO>=thres(2);
                                FRVO_control_movTW.low(itwin,ncells_ana)  = nanmean(neurFR_EO(ind_low));
                                FRVO_control_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR_EO(ind_mid));
                                FRVO_control_movTW.high(itwin,ncells_ana) = nanmean(neurFR_EO(ind_high));   
                            elseif SOslope_icell == -1 %CJB
                                thres = quantile(OfferValueA_EO,2);
                                ind_low  = OfferValueA_EO< thres(1);
                                ind_mid  = OfferValueA_EO< thres(2) & OfferValueA_EO>=thres(1);
                                ind_high = OfferValueA_EO>=thres(2);
                                FRVO_control_movTW.low(itwin,ncells_ana)  = nanmean(neurFR_EO(ind_low));
                                FRVO_control_movTW.mid(itwin,ncells_ana)  = nanmean(neurFR_EO(ind_mid));
                                FRVO_control_movTW.high(itwin,ncells_ana) = nanmean(neurFR_EO(ind_high)); 
                            end      
                        end % if JC or SO   
                        %
                        allneuroFR_OE_icell(:,itwin) = neurFR_OE; % each row: each trial; each colomn: each time window                    
                    end %for itwin        

                                        %
                    if JCslope_icell == 1 % CJA
                        valuerange_all(:,ncells_ana) = rho_all(:,ncells_ana).*(max(abs(neuract_all(:,2))) - min(abs(neuract_all(:,2))));    
                    elseif JCslope_icell == -1 %CJB
                        valuerange_all(:,ncells_ana) = (max(abs(neuract_all(:,3))) - min(abs(neuract_all(:,3)))); 
                    end
                    
                    
                    % %
                    % average firing rate across critical time window, rather than using moving time windows
                    if isequal(alignment,'off')
                    % average firing rate across critical time window, and do the correlation 
                    if isequal(JCorSO,'JC')
                        critTW = [250,350];
                        ind_TW = timepoints>=critTW(1) & timepoints<=critTW(2);
                        averCritFR = nanmean(allneuroFR_OE_icell(:,ind_TW),2); % each row: each trial; one colomn: average across critical TW
                        % calculate correlation between V(O) and average FR across critical TW
                        if JCslope_icell == 1 % CJA
                            [rr,pp] = corr(OfferValueB,averCritFR);
                            aa = polyfit(OfferValueB,averCritFR,1);
                            a0 = aa(2);
                            a1 = aa(1);
                        elseif JCslope_icell == -1 %CJB
                            [rr,pp] = corr(OfferValueA,averCritFR);
                            aa = polyfit(OfferValueA,averCritFR,1);
                            a0 = aa(2);
                            a1 = aa(1);
                        end
                        averCritFR_vs_VO_RR(:,ncells_ana) = rr;
                        averCritFR_vs_VO_PP(:,ncells_ana) = pp;
                        averCritFR_vs_VO_a0(:,ncells_ana) = a0;
                        averCritFR_vs_VO_a1(:,ncells_ana) = a1;
                        corr_VOandVE(:,ncells_ana) = corr(OfferValueA,OfferValueB);
                    elseif isequal(JCorSO,'SO')
                        % critTW = [1150,1450]; % after offer2 onset
                        critTW = [750,1050];  % before offer2 onset
                        % critTW = [750,950];  % before offer2 onset
                        ind_TW = timepoints>=critTW(1) & timepoints<=critTW(2);
                        averCritFR = nanmean(allneuroFR_OE_icell(:,ind_TW),2); % each row: each trial; one colomn: average across critical TW
                        % calculate correlation between V(O) and average FR across critical TW
                        if SOslope_icell == 1 % CJA
                            [rr,pp] = corr(OfferValueB_OE,averCritFR);
                            aa = polyfit(OfferValueB_OE,averCritFR,1);
                            a0 = aa(2);
                            a1 = aa(1);
                        elseif SOslope_icell == -1 %CJB
                            [rr,pp] = corr(OfferValueA_OE,averCritFR);
                            aa = polyfit(OfferValueA_OE,averCritFR,1);
                            a0 = aa(2);
                            a1 = aa(1);
                        end
                        averCritFR_vs_VO_RR(:,ncells_ana) = rr;
                        averCritFR_vs_VO_PP(:,ncells_ana) = pp;
                        averCritFR_vs_VO_a0(:,ncells_ana) = a0;
                        averCritFR_vs_VO_a1(:,ncells_ana) = a1;
                        corr_VOandVE(:,ncells_ana) = corr(OfferValueA_OE,OfferValueB_OE);
                    end
                    end                      
                end %for icell 
                
                % save results
                if ~JCSOclassifySep
                    filename = ['pop_circuit_inhibition_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,dobhvlogit];
%                     filename = ['pop_circuit_inhibition_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,'_neworderbias'];
                elseif JCSOclassifySep
                    filename = ['pop_circuit_inhibition_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,'_SOclassonly',dobhvlogit];
%                     filename = ['pop_circuit_inhibition_bothJCandSO_TT_',monkey_ana,'_',cellclassnames{iclass},'_',JCorSO,'_',alignment,'_SOclassonly_neworderbias'];      
                end
                eval(['save ',filename ' monkeyname_list cellnames_list tuningslopes_list ncells_ana rho_all orderbias_all Steepness_all valuerange_all timepoints '...
                       ' FRvsVO_RR_movTW FRvsVO_PP_movTW FRvsVO_a0_movTW FRvsVO_a1_movTW FRVO_movTW FRVE_movTW FRVO_control_movTW FRVE_control_movTW'...
                       ' a1vsOrderbias_RR_movTW PPvsSteep_PP_movTW averCritFR_vs_VO_RR averCritFR_vs_VO_PP averCritFR_vs_VO_a0 averCritFR_vs_VO_a1 corr_VOandVE ' ])
                
            end % try catch
            
            if isequal(monkey_plot,'both')
                ind_monplot = logical(ones(size(monkeyname_list)));
            elseif isequal(monkey_plot,'Gervinho')
                ind_monplot = ismember(monkeyname_list,'G');
            elseif isequal(monkey_plot,'Juan')
                ind_monplot = ismember(monkeyname_list,'J');    
            end
            
            % PLOT - critical time window
            if 0
            % calculate correlation between FRvsVO relationship and Steepness - critical time window
                ind_nan = isnan(Steepness_all);
                ind_nan = ind_nan | isnan(averCritFR_vs_VO_RR);
                ind_good = ~ind_nan & Steepness_all<15 & ind_monplot;
                %
                figure;set(gcf,'position',[110 65 850 850], 'PaperPositionMode','auto')
                subplot(2,2,1)
                XXX = (Steepness_all(:,ind_good))';
                YYY = averCritFR_vs_VO_RR(ind_good)';
                hold on; plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
                aa = polyfit(XXX,YYY,1);
                XX = [0.5, 8];
                YY = aa(1)*XX+aa(2);
                hold on; plot(XX,YY,'b-','LineWidth',3);
                Sigma_ell = cov(XXX, YYY);
                mu_ell(1) = nanmean(XXX);
                mu_ell(2) = nanmean(YYY);       
                hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                title('correlation between correlation of critical FR vs V(O) and steepness');
                xlabel('steepness');ylabel('corr');
                [RR_Pea, PP_Pea] = corr(XXX, YYY, 'type', 'Pearson');
                [RR_Spe, PP_Spe] = corr(XXX, YYY, 'type', 'Spearman');
                text(XX(1)+0.1,max(YYY)-0.1,{['corr(Pearson) = ',num2str(RR_Pea)]; ['p(Pearson) = ',num2str(PP_Pea)];...
                                             ['corr(Spearman) = ',num2str(RR_Spe)]; ['p(Spearman) = ',num2str(PP_Spe)];}, 'fontsize', 8);  
                %
                subplot(2,2,2)
                XXX = (Steepness_all(:,ind_good))';
                YYY = -log10(averCritFR_vs_VO_PP(ind_good))';
                hold on; plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
                plot(XXX,YYY,'ko','MarkerSize',8);
                aa = polyfit(XXX,YYY,1);
                XX = [0.5, 8];
                YY = aa(1)*XX+aa(2);
                hold on; plot(XX,YY,'b-','LineWidth',3);
                Sigma_ell = cov(XXX, YYY);
                mu_ell(1) = nanmean(XXX);
                mu_ell(2) = nanmean(YYY);       
                hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                title('correlation between -log10(p) of correlation of critical FR vs V(O) and steepness');
                xlabel('steepness');ylabel('-log10(p)');
                [RR_Pea, PP_Pea] = corr(XXX, YYY, 'type', 'Pearson');
                [RR_Spe, PP_Spe] = corr(XXX, YYY, 'type', 'Spearman');
                text(XX(1)+0.1,max(YYY)-0.1,{['corr(Pearson) = ',num2str(RR_Pea)]; ['p(Pearson) = ',num2str(PP_Pea)];...
                                             ['corr(Spearman) = ',num2str(RR_Spe)]; ['p(Spearman) = ',num2str(PP_Spe)];}, 'fontsize', 8);
                %
                subplot(2,2,3)
                XXX = (Steepness_all(:,ind_good))';
                YYY = averCritFR_vs_VO_a1(ind_good)';
                hold on; plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
                aa = polyfit(XXX,YYY,1);
                XX = [0.5, 8];
                YY = aa(1)*XX+aa(2);
                hold on; plot(XX,YY,'b-','LineWidth',3);
                Sigma_ell = cov(XXX, YYY);
                mu_ell(1) = nanmean(XXX);
                mu_ell(2) = nanmean(YYY);       
                hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                title('correlation between beta1(regression) of critical FR vs V(O) and steepness');
                xlabel('steepness');ylabel('beta1(regression)');
                [RR_Pea, PP_Pea] = corr(XXX, YYY, 'type', 'Pearson');
                [RR_Spe, PP_Spe] = corr(XXX, YYY, 'type', 'Spearman');
                text(XX(1)+0.1,max(YYY)-0.1,{['corr(Pearson) = ',num2str(RR_Pea)]; ['p(Pearson) = ',num2str(PP_Pea)];...
                                             ['corr(Spearman) = ',num2str(RR_Spe)]; ['p(Spearman) = ',num2str(PP_Spe)];}, 'fontsize', 8);
                %
                % average neurons within sessions
                steepness_sess = unique(Steepness_all','rows');
                averCritFR_vs_VO_PP_sess = [];
                averCritFR_vs_VO_RR_sess = [];
                averCritFR_vs_VO_a1_sess = [];
                nsess = length(steepness_sess);
                for isess = 1:nsess
                    steepness_isess = steepness_sess(isess);
                    ind_isess = steepness_sess == steepness_isess;
                    averCritFR_vs_VO_PP_sess(isess,:) = nanmean(averCritFR_vs_VO_PP(ind_isess));
                    averCritFR_vs_VO_RR_sess(isess,:) = nanmean(averCritFR_vs_VO_RR(ind_isess));
                    averCritFR_vs_VO_a1_sess(isess,:) = nanmean(averCritFR_vs_VO_a1(ind_isess));
                end
                ind_nan = isnan(steepness_sess);
                ind_nan = ind_nan | isnan(averCritFR_vs_VO_RR_sess) | isnan(averCritFR_vs_VO_PP_sess);
                ind_good = ~ind_nan;
                %
                figure;set(gcf,'position',[110 65 850 850], 'PaperPositionMode','auto')
                subplot(2,2,1)
                XXX = (steepness_sess(ind_good));
                YYY = averCritFR_vs_VO_RR_sess(ind_good);
                hold on; plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
                aa = polyfit(XXX,YYY,1);
                XX = [0.5, 8];
                YY = aa(1)*XX+aa(2);
                hold on; plot(XX,YY,'b-','LineWidth',3);
                Sigma_ell = cov(XXX, YYY);
                mu_ell(1) = nanmean(XXX);
                mu_ell(2) = nanmean(YYY);       
                hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                title({['correlation between correlation of critical FR vs V(O) and steepenss'];['average with session']});
                xlabel('steepness');ylabel('corr');
                [RR_Pea, PP_Pea] = corr(XXX, YYY, 'type', 'Pearson');
                [RR_Spe, PP_Spe] = corr(XXX, YYY, 'type', 'Spearman');
                text(XX(1)+0.1,max(YYY)-0.1,{['corr(Pearson) = ',num2str(RR_Pea)]; ['p(Pearson) = ',num2str(PP_Pea)];...
                                             ['corr(Spearman) = ',num2str(RR_Spe)]; ['p(Spearman) = ',num2str(PP_Spe)];}, 'fontsize', 8);
                %
                subplot(2,2,2)
                XXX = (steepness_sess(ind_good));
                YYY = -log10(averCritFR_vs_VO_PP_sess(ind_good));
                hold on; plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
                plot(XXX,YYY,'ko','MarkerSize',8);
                aa = polyfit(XXX,YYY,1);
                XX = [0.5, 8];
                YY = aa(1)*XX+aa(2);
                hold on; plot(XX,YY,'b-','LineWidth',3);
                Sigma_ell = cov(XXX, YYY);
                mu_ell(1) = nanmean(XXX);
                mu_ell(2) = nanmean(YYY);       
                hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                title({['correlation between -log10(p) of correlation of critical FR vs V(O) and steepness'];['average with session']});
                xlabel('steepness');ylabel('-log10(p)');
                [RR_Pea, PP_Pea] = corr(XXX, YYY, 'type', 'Pearson');
                [RR_Spe, PP_Spe] = corr(XXX, YYY, 'type', 'Spearman');
                text(XX(1)+0.1,max(YYY)-0.1,{['corr(Pearson) = ',num2str(RR_Pea)]; ['p(Pearson) = ',num2str(PP_Pea)];...
                                             ['corr(Spearman) = ',num2str(RR_Spe)]; ['p(Spearman) = ',num2str(PP_Spe)];}, 'fontsize', 8);
                %
                subplot(2,2,3)
                XXX = (steepness_sess(ind_good));
                YYY = averCritFR_vs_VO_a1_sess(ind_good);
                hold on; plot(XXX,YYY,'ko','MarkerSize',8,'LineWidth',1);
                aa = polyfit(XXX,YYY,1);
                XX = [0.5, 8];
                YY = aa(1)*XX+aa(2);
                hold on; plot(XX,YY,'b-','LineWidth',3);
                Sigma_ell = cov(XXX, YYY);
                mu_ell(1) = nanmean(XXX);
                mu_ell(2) = nanmean(YYY);       
                hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                title({['correlation between beta1(regression) of critical FR vs V(O) and steepenss'];['average with session']});
                xlabel('steepness');ylabel('beta1(regression)');
                [RR_Pea, PP_Pea] = corr(XXX, YYY, 'type', 'Pearson');
                [RR_Spe, PP_Spe] = corr(XXX, YYY, 'type', 'Spearman');
                text(XX(1)+0.1,max(YYY)-0.1,{['corr(Pearson) = ',num2str(RR_Pea)]; ['p(Pearson) = ',num2str(PP_Pea)];...
                                             ['corr(Spearman) = ',num2str(RR_Spe)]; ['p(Spearman) = ',num2str(PP_Spe)];}, 'fontsize', 8);
            end
            
            % %
            % %
            if 1
            % calculate correlation between FRvsVO relationship and order bias of SO - critical time window
                if isequal(JCorSO,'SO')
                    %
                    % remove orderbias_all outliers
                    kout = 3;
                    mean_OB = nanmean(orderbias_all);
                    std_OB = nanstd(orderbias_all);
                    outlier_OB = [mean_OB-kout*std_OB, mean_OB+kout*std_OB];
%                     kout = 1.5;
%                     quantiles_OB = quantile(orderbias_all,[0.25 0.5 0.75]);
%                     IQR_OB = quantiles_OB(3) - quantiles_OB(1);
%                     outlier_OB = [quantiles_OB(1)-kout*IQR_OB, quantiles_OB(3)+kout*IQR_OB];
                    ind_goodorder = orderbias_all > outlier_OB(1)  & orderbias_all < outlier_OB(2);
                    %                   
                    ind_nan = isnan(orderbias_all);
                    ind_nan = ind_nan | isnan(averCritFR_vs_VO_RR);
                    ind_nanout = ~ind_nan;
                    %
                    % only neurons with significant correlation
                    ind_goodcorr = averCritFR_vs_VO_PP<0.05;
                    %
                    % ind_good = ind_nanout & ind_goodorder & ind_goodcorr & ind_monplot;
                    ind_good = ind_nanout & ind_goodorder & ind_monplot;
                    anacellnum =  sum([ind_good]); 
                    %
                    figure;set(gcf,'position',[110 65 850 850], 'PaperPositionMode','auto')
                    subplot(2,2,1)
                    XXX = (orderbias_all(:,ind_good))';
                    YYY = averCritFR_vs_VO_RR(ind_good)';
                    XX = [floor(min(XXX)),ceil(max(XXX))];
                    YY = [floor(min(YYY)),ceil(max(YYY))]; 
                    %
                    aa = polyfit(XXX,YYY,1);
                    Yfit = aa(1)*XX+aa(2);
                    %
                    hold on; plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
                    hold on; plot(XXX, YYY, 'ko','MarkerSize',8);
                    Sigma_ell = cov(XXX, YYY);
                    mu_ell(1) = nanmean(XXX);
                    mu_ell(2) = nanmean(YYY);       
                    hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                    % title('correlation between correlation of critical FR vs V(O) and order bias');
                    xlabel('order bias');ylabel('corr');
                    [RR_Pea, pp_Pea] = corr(XXX, YYY, 'type', 'Pearson');
                    [RR_Spe, pp_Spe] = corr(XXX, YYY, 'type', 'Spearman');
                    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
                    {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
                     ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
                     ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
                    box off
                    axis([XX YY])
                    axis square    
                    set(gca,'FontSize',14);  
                    %
                    
                    subplot(2,2,2)
                    XXX = (orderbias_all(:,ind_good))';
                    YYY = -log10(averCritFR_vs_VO_PP(ind_good))';                   
                    aa = polyfit(XXX,YYY,1);
                    XX = [floor(min(XXX)),ceil(max(XXX))];
                    YY = [floor(min(YYY)),ceil(max(YYY))]; 
                    %
                    aa = polyfit(XXX,YYY,1);
                    Yfit = aa(1)*XX+aa(2);
                    %
                    hold on; plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
                    hold on; plot(XXX, YYY, 'ko','MarkerSize',8);
                    Sigma_ell = cov(XXX, YYY);
                    mu_ell(1) = nanmean(XXX);
                    mu_ell(2) = nanmean(YYY);       
                    hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                    % title('correlation between -log10(p) of correlation of critical FR vs V(O) and order bias');
                    xlabel('order bias');ylabel('-log10(p)');
                    [RR_Pea, pp_Pea] = corr(XXX, YYY, 'type', 'Pearson');
                    [RR_Spe, pp_Spe] = corr(XXX, YYY, 'type', 'Spearman');
                    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
                    {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
                     ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
                     ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
                    box off
                    axis([XX YY])
                    axis square    
                    set(gca,'FontSize',14); 
                    
                    %
                    subplot(2,2,3)
                    XXX = (orderbias_all(:,ind_good))';
                    YYY = [averCritFR_vs_VO_a1(ind_good)]';   
                    % YYY = [averCritFR_vs_VO_a1(ind_good).*valuerange_all(ind_good)]';   
                    %
                    % remove YYY outlier
                    kout = 3;
                    mean_Y = nanmean(YYY);
                    std_Y = nanstd(YYY);
                    outlier_Y = [mean_Y-kout*std_Y, mean_Y+kout*std_Y];
%                     kout = 1.5;
%                     quantiles_Y = quantile(YYY,[0.25 0.5 0.75]);
%                     IQR_Y = quantiles_Y(3) - quantiles_Y(1);
%                     outlier_Y = [quantiles_Y(1)-kout*IQR_Y, quantiles_Y(3)+kout*IQR_Y];
                    ind_goodYYY = YYY > outlier_Y(1)  & YYY < outlier_Y(2);
                    anacellnum = sum(ind_goodYYY);
                    XXX = XXX(ind_goodYYY);
                    YYY = YYY(ind_goodYYY);
                    %
                    aa = polyfit(XXX,YYY,1);
                    XX = [floor(min(XXX)),ceil(max(XXX))];
                    YY = [floor(min(YYY)),ceil(max(YYY))]; 
                    %
                    aa = polyfit(XXX,YYY,1);
                    Yfit = aa(1)*XX+aa(2);
                    %
                    hold on; plot(XX,Yfit,'-','LineWidth',3, 'Color', [0.4 0.4 0.4]);
                    hold on; plot(XXX, YYY, 'ko','MarkerSize',8);
                    Sigma_ell = cov(XXX, YYY);
                    mu_ell(1) = nanmean(XXX);
                    mu_ell(2) = nanmean(YYY);       
                    hold on; plotErrorEllipse(mu_ell, Sigma_ell, 0.90);
                    % title('correlation between beta1(regression) of critical FR vs V(O) and order bias');
                    xlabel('order bias');
                    ylabel('\beta1(regression)');
                    % ylabel('\beta1 * \DeltaV');
                    [RR_Pea, pp_Pea] = corr(XXX, YYY, 'type', 'Pearson');
                    [RR_Spe, pp_Spe] = corr(XXX, YYY, 'type', 'Spearman');
                    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
                    {['Spearman: r=',num2str(RR_Spe,'%1.1g'), ', p=',num2str(pp_Spe,'%1.1g')]; ...
                     ['Pearson: r=',num2str(RR_Pea,'%1.1g'),', p=',num2str(pp_Pea,'%1.1g')];...
                     ['N = ',num2str(anacellnum),' cells']}, 'fontsize', 12);
                    box off
                    axis([XX YY])
                    axis square    
                    set(gca,'FontSize',14); 
                                             
                    %
                    subplot(2,2,4)
                    % XXX = averCritFR_vs_VO_a1(ind_good)';
                    XXX = [averCritFR_vs_VO_a1(ind_good).*valuerange_all(ind_good)]';   
                    % XXX = [averCritFR_vs_VO_a1(:).*valuerange_all(:)]';   
                    % XX = [-2, 2]; 
                    XX = [-10, 10]; 
                    YY = [0, 50];
                    % edges = [XX(1):0.1:XX(2)];
                    edges = [XX(1):1:XX(2)];
                    hold on; 
                    plot([0 0],YY, '--','Color', [0.5 0.5 0.5],'LineWidth',1); hold on;
                    histogram(XXX,edges,'FaceColor',[.4 .4 .4]); hold on;                   
                    [~, pp_ttest] = ttest(XXX);
                    [ pp_wil,  ~] = signrank(XXX);                    
                    % title('distrbution of beta1(regression) of critical FR vs V(O)');
                    % xlabel('\beta1(regression)');
                    xlabel('\beta1 * \DeltaV');
                    ylabel('cell #');
                    text(XX(1)+(XX(2)-XX(1))/15, YY(2)-(YY(2)-YY(1))/10,...
                        {['t test: p = ',num2str(pp_ttest,'%1.1g')];...
                         ['Wilcoxon: p = ',num2str(pp_wil,'%1.1g')];...
                         ['mean = ',num2str(mean(XXX),'%.2f')];...                                        
                         },'fontsize',12);
                    box off
                    axis([XX YY])
                    axis square    
                    set(gca,'FontSize',14); 
                    
                    % 
                    axes('position',[.02 .97 .2 .05]);
                    text(0,0,{['chosen juice cells'];['300ms before offer2 (-250ms~50ms)'];['Task 2 in Two Task']},'fontsize',10);
                    axis off
                end
            end
            
            
            % % PLOT - sliding time windows
            if 0
                % calculate correlation between FRvsVO relationship and order bias - sliding time windows
                for itwin = 1:ntwins
                    ind_nan = isnan(orderbias_all);
                    ind_nan = ind_nan | isnan(FRvsVO_a1_movTW(itwin,:));
                    ind_good = ~ind_nan & ind_monplot;
                    [a1vsOrderbias_RR_movTW(itwin,:), a1vsOrderbias_PP_movTW(itwin,:)] = corr(orderbias_all(:,ind_good)',FRvsVO_a1_movTW(itwin,ind_good)');
                end

                %
                figure;set(gcf,'position',[110 65 1150 550], 'PaperPositionMode','auto')
                subplot(1,4,1)
                hold on
                plot(timepoints, nanmean(FRVO_movTW.low,2),'-','LineWidth',3);  
                plot(timepoints, nanmean(FRVO_movTW.mid,2),'-','LineWidth',3);  
                plot(timepoints, nanmean(FRVO_movTW.high,2),'-','LineWidth',3); 
                legend({'low V(O)','middle V(O)', 'high V(O)'})
                if isequal(alignment,'off')
                if isequal(JCorSO,'JC')
                    plot([0 0], [4 8], 'k--');
                elseif isequal(JCorSO,'SO')
                    plot([0 0], [4 8], 'k--'); hold on
                    plot([500 500], [4 8], 'k--'); hold on
                    plot([1000 1000], [4 8], 'k--'); hold on
                    plot([1500 1500], [4 8], 'k--'); hold on
                end
                xlabel('time(ms)')
                ylabel('average firing rate(spec/s)')
                title({['Firing rate in three V(0) quantiles'];['CJ cells; OffOn Align']});
                end
                %
                subplot(1,4,2)
                plot(timepoints, nanmean(FRvsVO_RR_movTW,2),'-'); hold on
                if isequal(alignment,'off')
                if isequal(JCorSO,'JC')
                    plot([0 0], [-0.1 0.1], 'k--');
                elseif isequal(JCorSO,'SO')
                    plot([0 0], [-0.1 0.1], 'k--'); hold on
                    plot([500 500], [-0.1 0.1], 'k--'); hold on
                    plot([1000 1000], [-0.1 0.1], 'k--'); hold on
                    plot([1500 1500], [-0.1 0.1], 'k--'); hold on
                end
                xlabel('time(ms)')
                ylabel('corr of FR vs V(O)')
                title({['corr between FR and V(O)'];['CJ cells; OffOn Align in ',JCorSO, ' trials']});
                end
                %
                subplot(1,4,3)
                plot(timepoints, nanmean(-log10(FRvsVO_PP_movTW),2),'-'); hold on
                if isequal(alignment,'off')
                if isequal(JCorSO,'JC')
                    plot([0 0], [0 1], 'k--');
                elseif isequal(JCorSO,'SO')
                    plot([0 0], [0 1], 'k--'); hold on
                    plot([500 500], [0 1], 'k--'); hold on
                    plot([1000 1000], [0 1], 'k--'); hold on
                    plot([1500 1500], [0 1], 'k--'); hold on
                end
                xlabel('time(ms)')
                ylabel('-log10(p) of FR vs V(O)')
                title({['-log10(p) between FR and V(O)'];['CJ cells; OffOn Align in ',JCorSO, ' trials']});
                end
                %
                subplot(1,4,4)
                negRR = FRvsVO_RR_movTW<0;
                posRR = FRvsVO_RR_movTW>0;
                signiPP = FRvsVO_PP_movTW<0.05;
                signinegRR = double(negRR).*double(signiPP);
                signiposRR = double(posRR).*double(signiPP);
                plot(timepoints,sum(signiposRR,2)./size(signiposRR,2),'r-');
                hold on
                plot(timepoints,sum(signinegRR,2)./size(signinegRR,2),'g-');
                if isequal(alignment,'off')
                if isequal(JCorSO,'JC')
                    plot([0 0], [0 0.1], 'k--');
                elseif isequal(JCorSO,'SO')
                    plot([0 0], [0 0.1], 'k--'); hold on
                    plot([500 500], [0 0.1], 'k--'); hold on
                    plot([1000 1000], [0 0.1], 'k--'); hold on
                    plot([1500 1500], [0 0.1], 'k--'); hold on
                end
                xlabel('time(ms)')
                ylabel('cell # fractions')
                legend({['positive encoding'],['negative encoding']});
                end

                % % 
                figure;set(gcf,'position',[110 65 1150 550], 'PaperPositionMode','auto')
                subplot(1,2,1)
                plot(timepoints, a1vsOrderbias_RR_movTW,'-');  hold on
                if isequal(alignment,'off')
                if isequal(JCorSO,'JC')
                    plot([0 0], [-0.1 0.1], 'k--');
                elseif isequal(JCorSO,'SO')
                    plot([0 0], [-0.1 0.1], 'k--'); hold on
                    plot([500 500], [-0.1 0.1], 'k--'); hold on
                    plot([1000 1000], [-0.1 0.1], 'k--'); hold on
                    plot([1500 1500], [-0.1 0.1], 'k--'); hold on
                end
                xlabel('time(ms)')
                ylabel('corr of a1 (FRvsV(O)) vs order bias')
                title({['CJ cells; OffOn Align in ',JCorSO, ' trials']});
                end
                %
                subplot(1,2,2)
                plot(timepoints, -log(a1vsOrderbias_PP_movTW),'-');  hold on
                if isequal(alignment,'off')
                if isequal(JCorSO,'JC')
                    plot([0 0], [0 1], 'k--');
                elseif isequal(JCorSO,'SO')
                    plot([0 0], [0 1], 'k--'); hold on
                    plot([500 500], [0 1], 'k--'); hold on
                    plot([1000 1000], [0 1], 'k--'); hold on
                    plot([1500 1500], [0 1], 'k--'); hold on
                end
                xlabel('time(ms)')
                ylabel('-log10(p) of a1 (FRvsV(O)) vs order bias')
                title({['CJ cells; OffOn Align in ',JCorSO, ' trials']});
                end
            end
            
            
            % % PLOT - sliding time windows
            if 1
                figure;set(gcf,'position',[110 65 1150 450], 'PaperPositionMode','auto')
                subplot(1,1,1)
                hold on
                ind_tmp = timepoints>=-450 & timepoints<=2550;
                ind_good = ind_monplot;
                plot(timepoints(ind_tmp), nanmean(FRVE_movTW.low(ind_tmp,ind_good),2), '-','LineWidth',3,'Color',[0.0 0.0 0.4]);  
                plot(timepoints(ind_tmp), nanmean(FRVE_movTW.mid(ind_tmp,ind_good),2), '-','LineWidth',3,'Color',[0.0 0.2 0.8]);  
                plot(timepoints(ind_tmp), nanmean(FRVE_movTW.high(ind_tmp,ind_good),2),'-','LineWidth',3,'Color',[0.2 0.8 1.0]);                 
                plot(timepoints(ind_tmp), nanmean(FRVO_movTW.low(ind_tmp,ind_good),2), '-','LineWidth',3,'Color',[0.6 0.0 0.2]);  
                plot(timepoints(ind_tmp), nanmean(FRVO_movTW.mid(ind_tmp,ind_good),2), '-','LineWidth',3,'Color',[1.0 0.0 0.0]);  
                plot(timepoints(ind_tmp), nanmean(FRVO_movTW.high(ind_tmp,ind_good),2),'-','LineWidth',3,'Color',[1.0 0.6 0.6]);  
                %
                legend({ 'EO trials; V(E) = Q1', 'EO trials; V(E) = Q2', 'EO trials; V(E) = Q3' ...
                         'OE trials; V(O) = Q1', 'OE trials; V(O) = Q2', 'OE trials; V(O) = Q3' ...
                        })
                if isequal(alignment,'off')
                if isequal(JCorSO,'JC')
                    YY = [3 9];
                    plot([0 0], [4 9], 'k--');
                elseif isequal(JCorSO,'SO')
                    YY = [4.5,9];
                    eventpoints = [0, 500, 1000, 1500, 2000];
                    eventnames = {'offer1 on','offer1 off','offer2 on','offer2 off','wait off'};
                    nevents = length(eventpoints);
                    for ievent = 1:nevents
                        eventpoint = [eventpoints(ievent),eventpoints(ievent)];
                        plot(eventpoint, YY, 'k--','LineWidth', 0.5); 
%                         text(eventpoints(ievent),YY(1)+(sum(YY)/25),eventnames{ievent},'fontsize',14);
                    end
                    set(gca,'XTick',eventpoints,'XTickLabel',eventnames);
                    %
                    plot([750 1050],[YY(1) YY(1)],'k-','LineWidth',4);
                    text(700,YY(1)+0.25,{'target time window'},'fontsize',12);
                end
                % xlabel('time (ms)')
                ylabel('Firing rate (sp/s)')
                title({['Task 2 in Two Task']});
                axis([-450 2350 YY]);
                set(gca,'Fontsize',17)
                box off
                end
            end
            
            % 
            % % PLOT - sliding time windows control analysis
            if 1
                figure;set(gcf,'position',[110 65 1150 450], 'PaperPositionMode','auto')
                subplot(1,1,1)
                hold on
                ind_tmp = timepoints>=-450 & timepoints<=2550;
                ind_good = ind_monplot;
                plot(timepoints(ind_tmp), nanmean(FRVE_control_movTW.low(ind_tmp,ind_good),2), '-','LineWidth',3,'Color',[0.0 0.0 0.4]);  
                plot(timepoints(ind_tmp), nanmean(FRVE_control_movTW.mid(ind_tmp,ind_good),2), '-','LineWidth',3,'Color',[0.0 0.2 0.8]);  
                plot(timepoints(ind_tmp), nanmean(FRVE_control_movTW.high(ind_tmp,ind_good),2),'-','LineWidth',3,'Color',[0.2 0.8 1.0]);                 
                plot(timepoints(ind_tmp), nanmean(FRVO_control_movTW.low(ind_tmp,ind_good),2), '-','LineWidth',3,'Color',[0.6 0.0 0.2]);  
                plot(timepoints(ind_tmp), nanmean(FRVO_control_movTW.mid(ind_tmp,ind_good),2), '-','LineWidth',3,'Color',[1.0 0.0 0.0]);  
                plot(timepoints(ind_tmp), nanmean(FRVO_control_movTW.high(ind_tmp,ind_good),2),'-','LineWidth',3,'Color',[1.0 0.6 0.6]);  
                %
                legend({ 'OE trials; V(E) = Q1', 'OE trials; V(E) = Q2', 'OE trials; V(E) = Q3' ...
                         'EO trials; V(O) = Q1', 'EO trials; V(O) = Q2', 'EO trials; V(O) = Q3' ...
                        })
                if isequal(alignment,'off')
                if isequal(JCorSO,'JC')
                    YY = [3 9];
                    plot([0 0], [4 9], 'k--');
                elseif isequal(JCorSO,'SO')
                    YY = [4.5,9];
                    eventpoints = [0, 500, 1000, 1500, 2000];
                    eventnames = {'offer1 on','offer1 off','offer2 on','offer2 off','wait off'};
                    nevents = length(eventpoints);
                    for ievent = 1:nevents
                        eventpoint = [eventpoints(ievent),eventpoints(ievent)];
                        plot(eventpoint, YY, 'k--','LineWidth', 0.5); 
%                         text(eventpoints(ievent),YY(1)+(sum(YY)/25),eventnames{ievent},'fontsize',14);
                    end
                    set(gca,'XTick',eventpoints,'XTickLabel',eventnames);
                    %
                    plot([750 1050],[YY(1) YY(1)],'k-','LineWidth',4);
                    text(700,YY(1)+0.25,{'target time window'},'fontsize',12);
                end
                % xlabel('time (ms)')
                ylabel('Firing rate (sp/s)')
                title({['Task 2 in Two Task']});
                axis([-450 2350 YY]);
                set(gca,'Fontsize',17)
                box off
                end
            end
            
        end %for ialign
%         %
%         eval(['Steepness_all_sumJCSO.',JCorSO,' = steepness_sess'';'])
%         %
%         eval(['averCritFR_vs_VO_RR_sumJCSO.',JCorSO,' = averCritFR_vs_VO_RR_sess'';'])
%         eval(['averCritFR_vs_VO_PP_sumJCSO.',JCorSO,' = averCritFR_vs_VO_PP_sess'';'])    
        %
        eval(['Steepness_all_sumJCSO.',JCorSO,' = Steepness_all;'])
        eval(['orderbias_all_sumJCSO.SO = orderbias_all;'])
        %
        eval(['averCritFR_vs_VO_RR_sumJCSO.',JCorSO,' = averCritFR_vs_VO_RR;'])
        eval(['averCritFR_vs_VO_PP_sumJCSO.',JCorSO,' = averCritFR_vs_VO_PP;']) 
        eval(['averCritFR_vs_VO_a1_sumJCSO.',JCorSO,' = averCritFR_vs_VO_a1;'])
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
