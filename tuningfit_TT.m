
function [datafit] = tuningfit_TT(cellname, anovastats, anovastats_both, atleast_nntrials, modelflag, differentRho) %#ok<STOUT,INUSL>
%
% Substitutes the old function do_tuningfit. 
%
% This function makes the linear regressions of neuronal responses on each
% of the many models defined in get_models.m the results are reported for
% different confidence intervals of nonzero regression slope. The parameter
% worththefit_thresholds is there to avoid regressing responses that do not
% pass the anova threshold criterion. 
%
% with modelflag = 'complete', all models are analyzed. 
% with modelflag = 'selected', only selected models are analyzed.

%
% author: camillo, september 2004
% revisions:	january 2005
% 				september 2005
% 				march 2007 (added start point to fit)
% 				september 2009
%				september 2016 for SO
%               october 2018 for TT - WS

if 0
	clear all %#ok<*UNRCH>
	tic
	cellname = 'L091207d11';
% 	modelflag = 'complete';
	modelflag = 'selected';
	atleast_nntrials = 3;			% min num of trials per trial type, for tuningfit
	anovastats = anovastats_TT(cellname, atleast_nntrials);
    anovastats_both = anovastats_both_TT(cellname, atleast_nntrials);
end
fitmode = 'tuningfit';

disp(['   ... fitting tuning of cell ',cellname,' (',fitmode,', ',modelflag,')'])

session = cellname(1:8); readsession_TT %#ok<NASGU>

filename = [dirroot, cellname, '_bhvParams']; eval(['load ',filename])

JCSO = {'SO','JC'};

for iJCorSO = 1:length(JCSO)
    JCorSO = JCSO{iJCorSO};
    
    [timewindows, ~] = get_timewindows_TT(JCorSO);
    % [timewindows ~] = get_timewindows(cellname);
    ntwins = length(timewindows);

    %selected models only?
    if 		isequal(modelflag,'complete')
        Rsqtouse = 'Rsq';
    elseif	isequal(modelflag,'selected')
        if isequal(JCorSO,'JC')
            selectedmodels = {'nA_off','nB_off','chosenvalue','Ach'};	%selected models
        elseif isequal(JCorSO,'SO')
            selectedmodels = {'offerAvalue1','offerBvalue1','chosenvalue1','Ach'};	%selected models
        end
        Rsqtouse = 'Rsq';
        cellclass = [];
    else
        disp('what models should be analyzed?'); dummy
    end

    %settings
    nonzero_confint = [ .90 .95 .99];	% nonzero slope confidence interval
    worththefit_threshold = 0.05;  	% threshold on anova to execute the fit % or 0.001 or 0.01 or 0.05  % change to 0.05 on 08/19/2019

    verbose = 0;	% verbose = 1 for one figure plot
    %				% verbose = 2 for separate plots


    %load tuning, sigmoidfit
    filename = [dirroot, cellname, '_tuning'];	eval(['load ',filename])
    try		filename = [dirroot, cellname, '_psyphycell']; eval(['load ',filename])
    catch,	psyphycell.sigmoidfit = sigmoidfit_TT(cellname); %#ok<*CTCH>
    end
    if isequal(JCorSO,'JC')
        pairnames = {'AB'};
        npairs = 1;
    elseif isequal(JCorSO,'SO')
        % pairnames = {'AB','BA','ABA'}; %fieldnames(tuning);
        pairnames = {'ABA'}; %fieldnames(tuning);
        npairs = 1;
    end

    %
    %keyboard
    for ipair = npairs:-1:1
    % 	pair = pairs(ipair);
        pairname = pairnames{ipair};
        
        if isequal(JCorSO,'JC') 
            eval(['table01 = tuning.',JCorSO,'.',pairname,'.table01;']);
        elseif isequal(JCorSO,'SO') 
            eval(['table01 = tuning.',JCorSO,'.',pairname,'.table02;']);
        end
            
		%eval(['table01 = tuning.',JCorSO,'.',pairname,'.table01;']);
        eval(['anovatab = anovastats.',JCorSO,'.',pairname,'.bytrialtype;'])
        eval(['anovatab_both = anovastats_both.JCSO.bytrialtype;']);
        %keyboard
        %get the models
        if differentRho
            eval(['models_all = get_models_',JCorSO,'(pairname, table01, psyphycell.sigmoidfit.',JCorSO,');']);
        elseif ~differentRho
            psyphycell.sigmoidfit.JC{3}(1) = (psyphycell.sigmoidfit.JC{3}(1)+psyphycell.sigmoidfit.SO{3}(1))/2;
            psyphycell.sigmoidfit.SO{3}(1) = (psyphycell.sigmoidfit.JC{3}(1)+psyphycell.sigmoidfit.SO{3}(1))/2;
            eval(['models_all = get_models_',JCorSO,'(pairname, table01, psyphycell.sigmoidfit.',JCorSO,');']);
        end
            
%       keyboard
        if 		isequal(modelflag,'complete')
            models = [];
            models = models_all;
        elseif	isequal(modelflag,'selected')
            models = [];
            for imodel = 1:size(selectedmodels,2)
                eval(['models.',selectedmodels{imodel},' = models_all.',selectedmodels{imodel},';'])
            end
        end
        modelnames = fieldnames(models);
        nmodels = size(modelnames,1);
	
        %initialize
        Rsq = NaN*ones(nmodels, length(timewindows));
        Rsqadj = NaN*ones(nmodels, length(timewindows));
        for iint = 1:length(nonzero_confint)
            nonzero_all{iint} = NaN*ones(nmodels, length(timewindows)); 
        end
        intercept	= NaN*ones(nmodels, length(timewindows));
        slope		= NaN*ones(nmodels, length(timewindows));
	
%       for iwin = 2
        for iwin = 1:length(timewindows)
            twin = timewindows{iwin}{1};
%             if isequal(twin,'postoffer')
%                 pval_current = min(anovatab.pval(iwin,1),anovatab_both.pval(1,1));
%             elseif isequal(twin,'postoffer1')
%                 pval_current = min(anovatab.pval(iwin,1),anovatab_both.pval(1,1));
%             elseif isequal(twin,'postoffer2')
%                 pval_current = min(anovatab.pval(iwin,1),anovatab_both.pval(1,1));
%             elseif isequal(twin,'postjuice') 
%                 pval_current = min(anovatab.pval(iwin,1),anovatab_both.pval(2,1));
%             else
%                 pval_current = anovatab.pval(iwin,1); 
%             end
            
            if isequal(JCorSO,'JC') 
                pval_threshold = anovatab.pval(iwin,1)<worththefit_threshold | iwin==2 | iwin==3 | iwin==6 | iwin==7 ; 
            elseif isequal(JCorSO,'SO') 
                pval_threshold = anovatab.pval(iwin,1)<worththefit_threshold | iwin==2 | iwin==4 | iwin==8 ; 
            end

            % if pval_current<worththefit_threshold
            if pval_threshold
            %	keyboard	
                %fit the models
                for imodel = 1:nmodels
                    [trialtypes, yy, ind_trialtype] = get_neuract_forfit_TT(tuning, pairname, twin, atleast_nntrials, JCorSO);
	%               keyboard
                    % TEST TO REMOVE FORCED CHOICES % SB
%                   trialtypes=trialtypes(~sum(trialtypes==0,2),:);
%                   ind_trialtype=logical(ind_trialtype(~sum(trialtypes==0,2)));
%                   yy=yy(~sum(trialtypes==0,2));
                    yy = yy(:,1);
                    modelname = modelnames{imodel};
                    xx_model = get_xxmodel_forfit(models, modelname, trialtypes);
                    xx = xx_model(ind_trialtype,end);
		
				
                    try
                        %linear fit
                        warning off
                        g = fittype('m*x+q','coeff',{'m','q'});
                        [linfit,goodness] = fit(xx,yy,g);
                        intercept(imodel,iwin) = linfit.q;
                        slope(imodel,iwin) = linfit.m;
                        for iint = 1:length(nonzero_confint)
                            interval = confint(linfit,nonzero_confint(iint));
                            if prod(sign(interval(:,1)))==1, nonzero = 1; else nonzero = 0; end
                            nonzero_all{iint}(imodel,iwin) = nonzero; 
                        end
                        Rsq(imodel,iwin) = goodness.rsquare;
                        Rsqadj(imodel,iwin) = goodness.adjrsquare;
	
                        if verbose, plot_tuningfit; end
				
                    catch
                        keyboard
                    end
				
                end %for imodel
			
            end %if worththefit
		
        end %for iwin
	
        %if selectedmodels, classifiy responses
        if 	isequal(modelflag, 'selected')
            clear cellclass
            cellclass.specs.fitmode = fitmode;
            cellclass.specs.Rsqtouse = Rsqtouse;
            cellclass.specs.anovathresh = worththefit_threshold;
            for iint = 1:length(nonzero_confint)
                nonzero_aux = nonzero_all{iint};		nonzero_aux(isnan(nonzero_aux)) = -1;
                eval(['Rsq_aux = ',Rsqtouse,';']);		Rsq_aux(isnan(Rsq_aux)) = -1; %#ok<AGROW>
                aux = nonzero_aux.*nonzero_aux.*Rsq_aux;			% nonzero_aux is multiplied 2x so that NaN -> -1 
                celclas = zeros(1,size(aux,2));						% cellclass = 0 by default
			%   
                [val,row] = max(aux,[],1);
                ind = val>0;	celclas(ind) = row(ind);			% cellclass
                ind = val==-1;	celclas(ind) = -1;	%#ok<NASGU>		% cellclass = -1 if NaN
			%
                pstr = ['p',num2str(100*nonzero_confint(iint))];
                eval(['cellclass.',pstr,' = celclas;'])
            end
        end
	
        %output
        for imodel = 1:nmodels
            modelname = modelnames{imodel};
            eval(['modellabs{imodel,1} = models.',modelname,'.modellab;'])
            eval(['modelclass{imodel,1} = models.',modelname,'.modelclass;'])
        end
	
        eval(['datafit.',JCorSO,'.',pairname,'.modellabs = modellabs;'])
        eval(['datafit.',JCorSO,'.',pairname,'.modelclass = modelclass;'])
        eval(['datafit.',JCorSO,'.',pairname,'.Rsq = Rsq;'])
        eval(['datafit.',JCorSO,'.',pairname,'.Rsqadj = Rsqadj;'])
        for iint = 1:length(nonzero_confint)
            pstr = ['p',num2str(100*nonzero_confint(iint))];
            eval(['datafit.',JCorSO,'.',pairname,'.nonzero.',pstr,' = nonzero_all{iint};'])
        end
        eval(['datafit.',JCorSO,'.',pairname,'.intercept = intercept;'])
        eval(['datafit.',JCorSO,'.',pairname,'.slope = slope;'])
	
        if 	isequal(modelflag, 'selected')
            eval(['datafit.',JCorSO,'.',pairname,'.cellclass = cellclass;'])
        end	
	
    end %for ipair
end %for iJCorSO