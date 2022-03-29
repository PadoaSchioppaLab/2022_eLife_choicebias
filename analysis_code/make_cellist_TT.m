% make_cellist_TT
% This code used to summarize the tuning property and cell classes for all
% OFC neurons

clearvars

%hyper parameters
monkey_ana = 'both'; % 'Gervinho', 'Juan', 'both'
eval(['sessionlist_',monkey_ana]);

brainarea = 'OFC';
pairname_JC = 'AB';
pairname_SO = 'ABA';

anovathresh = .01; % 0.001 or 0.01 or 0.05
anovaTWs_JC = [2,3,6,7]; % or [1:7]
anovaTWs_SO = [2,4,8];   % or [1:9]

impose_samesign = 1;

mintrialnum =  100; % 200, 160, 125, 100

% setup name for saving
filesavename = ['TTcellist_OFC_both'];

%% for JC trials
tic
JCcellist = [];

modelclasses = [12 13 6 14]; 
[timewindows,~] = get_timewindows_TT('JC');
ntwins = size(timewindows,2);
selected_twins = [2 3 6 7];

%initialize
ncells = 0;
infocell = [];
fr = []; 
cv = [];

for isession =	1:size(sessions,1)
	session = sessions{isession};
	readsession_TT
	
	%cells
	for icell = 1:size(cells_td,1)
		cellname = [num2str(cells_td(icell,1)),num2str(cells_td(icell,2))];
		
        try
            if checksession([session,cellname],.8), continue, end
        catch
            continue
        end
        
        %%%%%%%%%%%%
        filename = [dirroot,session,cellname,'_data'];
        eval(['load ',filename])
        clear celldata celldataerror psyphydata trace trialRecord
        if size(goodTrials_JC,1)<mintrialnum | size(goodTrials_SO,1)<mintrialnum % 200 or 160 or 125 or 100
            continue
        end
        %%%%%%%%%%%%
        
        
        if ~isequal(arearead_TT(session, cellname), brainarea), continue, end
		ncells = ncells+1;
		if ~rem(ncells,50), disp(['    processing cell number ',num2str(ncells)]); end
		
        
        filename = [dirroot,session,cellname,'_cellstats'];
        eval(['load ',filename])
               
        
		%find class across time windows
		eval(['ind_JC = cellstats.anovastats.JC.',pairname_JC,'.bytrialtype.pval(anovaTWs_JC,1)<anovathresh;'])
        eval(['ind_SO = cellstats.anovastats.SO.',pairname_SO,'.bytrialtype.pval(anovaTWs_SO,1)<anovathresh;'])
        infocell(ncells).passANOVA	= ~isempty(find(ind_JC,1))| ~isempty(find(ind_SO,1));
		if isempty(find(ind_JC,1)) && isempty(find(ind_SO,1)) % ANOVA threshold on either JC or SO 
        % if isempty(find(ind_JC,1)) % ANOVA threshold on JC only 
			classified = 0; 
			cellclass = -1; % do not pass ANOVA
			slopesign = 0;
        else
            eval(['Rsq = cellstats.tuningfit.JC.',pairname_JC,'.Rsq(modelclasses,:);'])
			eval(['nonzero = cellstats.tuningfit.JC.',pairname_JC,'.nonzero.p95(modelclasses,:);'])
			eval(['intercept = cellstats.tuningfit.JC.',pairname_JC,'.intercept(modelclasses,:);'])
			eval(['slope = cellstats.tuningfit.JC.',pairname_JC,'.slope(modelclasses,:);'])
            
            %
			if ~impose_samesign		%we dont impose the same slope sign
				Rsq = Rsq.*nonzero;
                Rsq_all = Rsq(:,selected_twins)';
				[classified, clas] = max(sum(Rsq_all,1));
                if ismember(clas,0)
                    infocell(ncells).Rsq_all = 0;
                else
                    infocell(ncells).Rsq_all = Rsq_all(:,clas)';
                end
			else					%we impose the same slope sign
				Rsq_sign = [(Rsq.*nonzero.*(slope>0));(Rsq.*nonzero.*(slope<0))];
                Rsq_all = Rsq_sign(:,selected_twins)';
                Rsq_all(isnan(Rsq_all))=0;
				[classified, clas] = max(sum(Rsq_all,1));
                if ismember(clas,0)
                    infocell(ncells).Rsq_all = 0;
                else
                    infocell(ncells).Rsq_all = Rsq_all(:,clas)';
                end
				if clas>4
					clas = clas-4;
					slopesign = -1;
				else
					slopesign = +1;
				end
			end
			%
			%impose 1,2,3 convention
			if ~classified, cellclass = 0; slopesign = 0;  %untuned
			elseif	ismember(clas,[1 2]),	cellclass = 1;	%offer value
			elseif	ismember(clas,3),		cellclass = 2;	%chosen value
			elseif	ismember(clas,4),		cellclass = 3;	%taste
			end
			%			
		end %if ~isempty(ind_JC)
		
		% add valrange info
        beta = cellstats.psyphycell.sigmoidfit.JC{2};
		logRho = -beta(1)/beta(2);
		relvalue = exp(logRho);
		%
		filename = [dirroot,session,cellname,'_tuning'];
		eval(['load ',filename])
		eval(['act = tuning.JC.',pairname_JC,'.neuract.bytrial.preoffer;'])
		preoffer_fr = mean(act(:,end));
		%
		%we express slopes and valranges in uB
		if classified
			[valrange, valminmax] = get_valrange_cell(parsession, clas, tuning, relvalue,'JC');
			%intercepts, slopes and actranges for the sel_twins
			non0 = nonzero(clas,selected_twins);
			non0(non0==0) = nan;
			intercepts = intercept(clas,selected_twins) .* non0;
			slopes = slope(clas,selected_twins) .* non0;
			actranges = slopes*valrange;
			%
			%if B preferred, do all the necessary adjustments
			if relvalue<1
				relvaluememo = 1;
				relvalue = 1./relvalue;
				if ismember(clas,[1 2])
					clas = setdiff([1,2],clas);
				end
				valminmax = valminmax*relvalue;
				valrange = valrange*relvalue;
				slopes = slopes/relvalue;
			end
		else
			valminmax = nan;
			valrange = nan;
			intercepts = nan;
			slopes = nan;
			actranges = nan;
		end
		unival = logRho;
		
		if ~classified, clas = 0; end
		
		infocell(ncells).cellname	= [session, cellname];
		infocell(ncells).exp		= 'JC';
		infocell(ncells).cellclass	= cellclass;
		infocell(ncells).subclass	= clas;			%adjusted for cases with B preferred
		infocell(ncells).slopesign	= slopesign;
        
		%
		infocell(ncells).valueinfo.univalue		= unival;
		infocell(ncells).valueinfo.preoffer_fr	= preoffer_fr;
		infocell(ncells).valueinfo.valminmax	= valminmax;
		infocell(ncells).valueinfo.valrange		= valrange;
		infocell(ncells).valueinfo.intercept	= intercepts;
		infocell(ncells).valueinfo.slope		= slopes;
		infocell(ncells).valueinfo.actrange		= actranges;
		
		%load isi info
		filename = [dirroot,session,cellname,'_tuning'];
		eval(['load ',filename])
		fr = nan(1,size(selected_twins,2));
		for ii = 1:size(selected_twins,2)
			act = [];
			iwin = selected_twins(ii);
			twin = timewindows{iwin}{1};
			eval(['act = [act; tuning.JC.',pairname_JC,'.neuract.bytrial.',twin,'];'])
			fr(ii) = mean(act(:,end));
		end
		infocell(ncells).fr = fr;
		
	end %for icell
end %for isession

%output
specs.ncells = ncells;
twins = {};
ii = 0;
for iwin = selected_twins
	ii = ii+1;
	twins(ii) = timewindows{iwin}(1);
end
specs.twins = twins;
specs.impose_samesign = impose_samesign;

JCcellist.infocell = infocell;
JCcellist.specs = specs;

toc

%
%% for SO trials
tic
SOcellist = [];

% % three time windows
modelclasses = [2 3 16 6 5 17 7 9 15 20 19]; 

[timewindows,~] = get_timewindows_TT('SO');
ntwins = size(timewindows,2);
selected_twins = [ 2 4 8 ]; % postoffer1 postoffer2 postjuice

%initialize
ncells = 0;
infocell = [];
fr = []; 
cv = [];
Rsq_all = [];

for isession =	1:size(sessions,1)
	session = sessions{isession};
	readsession_TT
	
	%cells
	for icell = 1:size(cells_td,1)
		cellname = [num2str(cells_td(icell,1)),num2str(cells_td(icell,2))];
		
        try
            if checksession([session,cellname],.8), continue, end
        catch
            continue
        end
        
        filename = [dirroot,session,cellname,'_data'];
        eval(['load ',filename])
        clear celldata celldataerror psyphydata trace trialRecord
        if size(goodTrials_JC,1)<mintrialnum | size(goodTrials_SO,1)<mintrialnum % 200 or 160 or 125 or 100
            continue
        end
                            
        if ~isequal(arearead_TT(session, cellname), brainarea), continue, end
		ncells = ncells+1;
		if ~rem(ncells,50), disp(['    processing cell number ',num2str(ncells)]); end
		
       
        filename = [dirroot,session,cellname,'_cellstats'];
        eval(['load ',filename])
        
            
		%find class across time windows
		eval(['ind_JC = cellstats.anovastats.JC.',pairname_JC,'.bytrialtype.pval(anovaTWs_JC,1)<anovathresh;'])
        eval(['ind_SO = cellstats.anovastats.SO.',pairname_SO,'.bytrialtype.pval(anovaTWs_SO,1)<anovathresh;'])
        infocell(ncells).passANOVA	= ~isempty(find(ind_SO,1)) | ~isempty(find(ind_JC,1));
		if isempty(find(ind_JC,1)) && isempty(find(ind_SO,1)) % ANOVA threshold on either JC or SO 
        % if isempty(find(ind_SO,1)) % ANOVA threshold on SO only 
			classified = 0; 
			cellclass = -1; % do not pass ANOVA
			slopesign = 0;
        else
            eval(['Rsq = cellstats.tuningfit.SO.',pairname_SO,'.Rsq(modelclasses,:);'])
			eval(['nonzero = cellstats.tuningfit.SO.',pairname_SO,'.nonzero.p95(modelclasses,:);'])
			eval(['intercept = cellstats.tuningfit.SO.',pairname_SO,'.intercept(modelclasses,:);'])
			eval(['slope = cellstats.tuningfit.SO.',pairname_SO,'.slope(modelclasses,:);'])
            
            %
			if ~impose_samesign		%we dont impose the same slope sign
				Rsq = Rsq.*nonzero;
                Rsq_all(:,1) = diag(Rsq(1:3,selected_twins)); % OVA
                Rsq_all(:,2) = diag(Rsq(4:6,selected_twins)); % OVB
                Rsq_all(:,3) = diag(Rsq(7:9,selected_twins)); % CV
                Rsq_all(:,4) = diag(Rsq([10 10 11],selected_twins)); % CJ
                
                [classified, clas] = max(sum(Rsq_all,1));
                if ismember(clas,0)
                    infocell(ncells).Rsq_all = 0;
                else
                    infocell(ncells).Rsq_all = Rsq_all(:,clas)';
                end
			else					%we impose the same slope sign
				Rsq_sign = [(Rsq.*nonzero.*(slope>0));(Rsq.*nonzero.*(slope<0))];
                Rsq_all(:,1) = diag(Rsq_sign(1:3,selected_twins)); % OVA+
                Rsq_all(:,2) = diag(Rsq_sign(4:6,selected_twins)); % OVB+
                Rsq_all(:,3) = diag(Rsq_sign(7:9,selected_twins)); % CV+
                Rsq_all(:,4) = diag(Rsq_sign([21 10 11],selected_twins)); % CJA
                Rsq_all(:,5) = diag(Rsq_sign(12:14,selected_twins)); % OVA-
                Rsq_all(:,6) = diag(Rsq_sign(15:17,selected_twins)); % OVB-
                Rsq_all(:,7) = diag(Rsq_sign(18:20,selected_twins)); % CV-
                Rsq_all(:,8) = diag(Rsq_sign([10 21 22],selected_twins)); % CJB
                
                Rsq_all(isnan(Rsq_all))=0;
				[classified, clas] = max(sum(Rsq_all,1));
                if ismember(clas,0)
                    infocell(ncells).Rsq_all = 0;
                else
                    infocell(ncells).Rsq_all = Rsq_all(:,clas)';
                end

				if clas>4
					clas = clas-4;
					slopesign = -1;
				else
					slopesign = +1;
				end
			end
			%
			%impose 1,2,3 convention
			if ~classified, cellclass = 0; slopesign = 0;  %untuned
			elseif	ismember(clas,[1 2]),	cellclass = 1;	%offer value
			elseif	ismember(clas,3),		cellclass = 2;	%chosen value
			elseif	ismember(clas,4),		cellclass = 3;	%taste
			end
			%			
		end %if ~isempty(ind)
		
		% add valrange info
        beta = cellstats.psyphycell.sigmoidfit.SO{2};
        % beta = cellstats.psyphycell.sigmoidfit.both{2};
		logRho = -beta(1)/beta(2);
		relvalue = exp(logRho);
		%
		filename = [dirroot,session,cellname,'_tuning'];
		eval(['load ',filename])
		eval(['act = tuning.SO.',pairname_SO,'.neuract.bytrial.preoffers;'])
		preoffer_fr = mean(act(:,end));
		%
		%we express slopes and valranges in uB
		if classified
			[valrange, valminmax] = get_valrange_cell(parsession, clas, tuning, relvalue, 'SO');
			%intercepts, slopes and actranges for the sel_twins
            if clas == 1     % OVA    
                clas_ind = [1:3];
            elseif clas == 2 % OVB    
                clas_ind = [4:6]; 
            elseif clas == 3 % CV    
                 clas_ind = [7:9]; 
            elseif clas == 4 % CJ    
                clas_ind = [10 10 11]; 
            end
			non0 = diag(nonzero(clas_ind,selected_twins));
			non0(non0==0) = nan;
			intercepts = diag(intercept(clas_ind,selected_twins)) .* non0;
			slopes = diag(slope(clas_ind,selected_twins)) .* non0;
			actranges = slopes*valrange;
			%
			%if B preferred, do all the necessary adjustments
			if relvalue<1
				relvaluememo = 1;
				relvalue = 1./relvalue;
				if ismember(clas,[1 2])
					clas = setdiff([1,2],clas);
				end
				valminmax = valminmax*relvalue;
				valrange = valrange*relvalue;
				slopes = slopes/relvalue;
			end
		else
			valminmax = nan;
			valrange = nan;
			intercepts = nan;
			slopes = nan;
			actranges = nan;
		end
		unival = logRho;
		
		if ~classified, clas = 0; end
		
		infocell(ncells).cellname	= [session, cellname];
		infocell(ncells).exp		= 'SO';
		infocell(ncells).cellclass	= cellclass;
		infocell(ncells).subclass	= clas;			%adjusted for cases with B preferred
		infocell(ncells).slopesign	= slopesign;
		%
		infocell(ncells).valueinfo.univalue		= unival;
		infocell(ncells).valueinfo.preoffer_fr	= preoffer_fr;
		infocell(ncells).valueinfo.valminmax	= valminmax;
		infocell(ncells).valueinfo.valrange		= valrange;
		infocell(ncells).valueinfo.intercept	= intercepts;
		infocell(ncells).valueinfo.slope		= slopes;
		infocell(ncells).valueinfo.actrange		= actranges;
		
		%load isi info
		filename = [dirroot,session,cellname,'_tuning'];
		eval(['load ',filename])
		fr = nan(1,size(selected_twins,2));
		for ii = 1:size(selected_twins,2)
			act = [];
			iwin = selected_twins(ii);
			twin = timewindows{iwin}{1};
			eval(['act = [act; tuning.SO.',pairname_SO,'.neuract.bytrial.',twin,'];'])
			fr(ii) = mean(act(:,end));
		end
		infocell(ncells).fr = fr;
		
	end %for icell
end %for isession

%output
specs.ncells = ncells;
twins = {};
ii = 0;
for iwin = selected_twins
	ii = ii+1;
	twins(ii) = timewindows{iwin}(1);
end
specs.twins = twins;
specs.impose_samesign = impose_samesign;

SOcellist.infocell = infocell;
SOcellist.specs = specs;

toc


%
%% for JC and SO trials together 
tic
JCSOcellist = [];

% % three time windows for both JC and SO
modelclasses_JC = [12 13 6 14];
modelclasses_SO = [2 3 16 6 5 17 7 9 15 20 19]; 

[timewindows_JC,~] = get_timewindows_TT('JC');
[timewindows_SO,~] = get_timewindows_TT('SO');

ntwins_JC = size(timewindows_JC,2);
ntwins_SO = size(timewindows_SO,2);
selected_twins_JC = [ 2 3 7 ]; % postoffer latedelay postjuice
selected_twins_SO = [ 2 4 8 ]; % postoffer1 postoffer2 postjuice

%initialize
ncells = 0;
infocell = [];
fr = []; 
cv = [];
Rsq_all_JC = [];
Rsq_all_SO = [];
Rsq_all = [];

for isession =	1:size(sessions,1)
	session = sessions{isession};
	readsession_TT
	
	%cells
	for icell = 1:size(cells_td,1)
		cellname = [num2str(cells_td(icell,1)),num2str(cells_td(icell,2))];
		
        try
            if checksession([session,cellname],.8), continue, end
        catch
            continue
        end
        
        filename = [dirroot,session,cellname,'_data'];
        eval(['load ',filename])
        clear celldata celldataerror psyphydata trace trialRecord
        if size(goodTrials_JC,1)<mintrialnum | size(goodTrials_SO,1)<mintrialnum % 200 or 160 or 125 or 100
            continue
        end
               
        
        if ~isequal(arearead_TT(session, cellname), brainarea), continue, end
		ncells = ncells+1;
		if ~rem(ncells,50), disp(['    processing cell number ',num2str(ncells)]); end
		
        
        filename = [dirroot,session,cellname,'_cellstats'];
        eval(['load ',filename])
 
		%find class across time windows
		eval(['ind_JC = cellstats.anovastats.JC.',pairname_JC,'.bytrialtype.pval(anovaTWs_JC,1)<anovathresh;'])
        eval(['ind_SO = cellstats.anovastats.SO.',pairname_SO,'.bytrialtype.pval(anovaTWs_SO,1)<anovathresh;'])
        infocell(ncells).passANOVA	= ~isempty(find(ind_JC,1)) | ~isempty(find(ind_SO,1));
		if isempty(find(ind_JC,1)) && isempty(find(ind_SO,1)) % ANOVA threshold on either JC or SO 
			classified = 0; % do not pass ANOVA
			cellclass = -1;
			slopesign = 0;
        else
            % JC
            eval(['Rsq_JC = cellstats.tuningfit.JC.',pairname_JC,'.Rsq(modelclasses_JC,:);'])
			eval(['nonzero_JC = cellstats.tuningfit.JC.',pairname_JC,'.nonzero.p95(modelclasses_JC,:);'])
			eval(['intercept_JC = cellstats.tuningfit.JC.',pairname_JC,'.intercept(modelclasses_JC,:);'])
			eval(['slope_JC = cellstats.tuningfit.JC.',pairname_JC,'.slope(modelclasses_JC,:);'])
            %SO
            eval(['Rsq_SO = cellstats.tuningfit.SO.',pairname_SO,'.Rsq(modelclasses_SO,:);'])
			eval(['nonzero_SO = cellstats.tuningfit.SO.',pairname_SO,'.nonzero.p95(modelclasses_SO,:);'])
			eval(['intercept_SO = cellstats.tuningfit.SO.',pairname_SO,'.intercept(modelclasses_SO,:);'])
			eval(['slope_SO = cellstats.tuningfit.SO.',pairname_SO,'.slope(modelclasses_SO,:);'])
            
            %
			if ~impose_samesign		%we dont impose the same slope sign
                % JC trials
                Rsq_JC = Rsq_JC.*nonzero_JC;
                Rsq_all_JC = Rsq_JC(:,selected_twins_JC)';
                % SO trials
                Rsq_SO = Rsq_SO.*nonzero_SO;                        
                Rsq_all_SO(:,1) = diag(Rsq_SO(1:3,selected_twins_SO)); % OVA
                Rsq_all_SO(:,2) = diag(Rsq_SO(4:6,selected_twins_SO)); % OVB
                Rsq_all_SO(:,3) = diag(Rsq_SO(7:9,selected_twins_SO)); % CV
                Rsq_all_SO(:,4) = diag(Rsq_SO([10 10 11],selected_twins_SO)); % CJ
                
                % consider both JC and SO
                Rsq_all = Rsq_all_JC + Rsq_all_SO;
				[classified, clas] = max(sum(Rsq_all,1));
			else					%we impose the same slope sign
				% JC trials
                Rsq_sign_JC = [(Rsq_JC.*nonzero_JC.*(slope_JC>0));(Rsq_JC.*nonzero_JC.*(slope_JC<0))];
                Rsq_all_JC = Rsq_sign_JC(:,selected_twins_JC)';
                Rsq_all_JC(isnan(Rsq_all_JC))=0;
                % SO trials
                Rsq_sign_SO = [(Rsq_SO.*nonzero_SO.*(slope_SO>0));(Rsq_SO.*nonzero_SO.*(slope_SO<0))];
                Rsq_all_SO(:,1) = diag(Rsq_sign_SO(1:3,selected_twins_SO)); % OVA+
                Rsq_all_SO(:,2) = diag(Rsq_sign_SO(4:6,selected_twins_SO)); % OVB+
                Rsq_all_SO(:,3) = diag(Rsq_sign_SO(7:9,selected_twins_SO)); % CV+
                Rsq_all_SO(:,4) = diag(Rsq_sign_SO([21 10 11],selected_twins_SO)); % CJA
                Rsq_all_SO(:,5) = diag(Rsq_sign_SO(12:14,selected_twins_SO)); % OVA-
                Rsq_all_SO(:,6) = diag(Rsq_sign_SO(15:17,selected_twins_SO)); % OVB-
                Rsq_all_SO(:,7) = diag(Rsq_sign_SO(18:20,selected_twins_SO)); % CV-
                Rsq_all_SO(:,8) = diag(Rsq_sign_SO([10 21 22],selected_twins_SO)); % CJB
                
                Rsq_all_SO(isnan(Rsq_all_SO))=0;
                % consider both JC and SO
                Rsq_all = Rsq_all_JC + Rsq_all_SO;
				[classified, clas] = max(sum(Rsq_all,1));
				if clas>4
					clas = clas-4;
					slopesign = -1;
				else
					slopesign = +1;
				end
			end
			%
			%impose 1,2,3 convention
			if ~classified, cellclass = 0; slopesign = 0;  %untuned
			elseif	ismember(clas,[1 2]),	cellclass = 1;	%offer value
			elseif	ismember(clas,3),		cellclass = 2;	%chosen value
			elseif	ismember(clas,4),		cellclass = 3;	%taste
			end
			
		end %if ~isempty(ind)
		
		% add valrange info
        % JC trials
        beta = cellstats.psyphycell.sigmoidfit.JC{2};
		logRho_JC = -beta(1)/beta(2);
		relvalue_JC = exp(logRho_JC);
        % SO trials
        beta = cellstats.psyphycell.sigmoidfit.SO{2};
		logRho_SO = -beta(1)/beta(2);
		relvalue_SO = exp(logRho_SO);
		%
		filename = [dirroot,session,cellname,'_tuning'];
		eval(['load ',filename])
        % JC trials
        eval(['act = tuning.JC.',pairname_JC,'.neuract.bytrial.preoffer;'])
		preoffer_fr_JC = mean(act(:,end));
        % SO trials
		eval(['act = tuning.SO.',pairname_SO,'.neuract.bytrial.preoffers;'])
		preoffer_fr_SO = mean(act(:,end));
        
        preoffer_fr = [preoffer_fr_JC;preoffer_fr_SO];
        %
		%we express slopes and valranges in uB
		if classified
            % JC trials
            [valrange_JC, valminmax_JC] = get_valrange_cell(parsession, clas, tuning, relvalue_JC,'JC');
			%intercepts, slopes and actranges for the sel_twins
			non0_JC = nonzero_JC(clas,selected_twins_JC);
			non0_JC(non0_JC==0) = nan;
			intercepts_JC = intercept_JC(clas,selected_twins_JC) .* non0_JC;
			slopes_JC = slope_JC(clas,selected_twins_JC) .* non0_JC;
			actranges_JC = slopes_JC*valrange_JC;
            %
            % SO trials
			[valrange_SO, valminmax_SO] = get_valrange_cell(parsession, clas, tuning, relvalue_SO, 'SO');
			%intercepts, slopes and actranges for the sel_twins
            if clas == 1     % OVA    
                clas_ind = [1:3]; 
            elseif clas == 2 % OVB    
                clas_ind = [4:6]; 
            elseif clas == 3 % CV    
                clas_ind = [7:9];
            elseif clas == 4 % CJ    
                clas_ind = [10 10 11]; 
            end
			non0_SO = diag(nonzero_SO(clas_ind,selected_twins_SO));
			non0_SO(non0_SO==0) = nan;
			intercepts_SO = diag(intercept_SO(clas_ind,selected_twins_SO)) .* non0_SO;
			slopes_SO = diag(slope_SO(clas_ind,selected_twins_SO)) .* non0_SO;
			actranges_SO = slopes_SO*valrange_SO;
			%
			%if B preferred, do all the necessary adjustments
			if (relvalue_JC + relvalue_SO)/2 < 1
				relvaluememo = 1;
				relvalue_JC = 1./relvalue_JC;
                relvalue_SO = 1./relvalue_SO;
				if ismember(clas,[1 2])
					clas = setdiff([1,2],clas);
				end
				valminmax_JC = valminmax_JC*relvalue_JC;
				valrange_JC = valrange_JC*relvalue_JC;
				slopes_JC = slopes_JC/relvalue_JC;
                valminmax_SO = valminmax_SO*relvalue_SO;
				valrange_SO = valrange_SO*relvalue_SO;
				slopes_SO = slopes_SO/relvalue_SO;
            end
            valminmax = [valminmax_JC';valminmax_SO'];
			valrange = [valrange_JC;valrange_SO];
			intercepts = [intercepts_JC;intercepts_SO'];
			slopes = [slopes_JC;slopes_SO'];
			actranges = [actranges_JC;actranges_SO'];
		else
			valminmax = [nan;nan];
			valrange = [nan;nan];
			intercepts = [nan;nan];
			slopes = [nan;nan];
			actranges = [nan;nan];
		end
		unival = [logRho_JC;logRho_SO];
		
		if ~classified, clas = 0; end
		
		infocell(ncells).cellname	= [session, cellname];
		infocell(ncells).exp		= 'JCSO';
		infocell(ncells).cellclass	= cellclass;
		infocell(ncells).subclass	= clas;			%adjusted for cases with B preferred
		infocell(ncells).slopesign	= slopesign;
		%
		infocell(ncells).valueinfo.univalue		= unival;
		infocell(ncells).valueinfo.preoffer_fr	= preoffer_fr;
		infocell(ncells).valueinfo.valminmax	= valminmax;
		infocell(ncells).valueinfo.valrange		= valrange;
		infocell(ncells).valueinfo.intercept	= intercepts;
		infocell(ncells).valueinfo.slope		= slopes;
		infocell(ncells).valueinfo.actrange		= actranges;
		
		%load isi info
		filename = [dirroot,session,cellname,'_tuning'];
		eval(['load ',filename])
		%JC trials
        fr_JC = nan(1,size(selected_twins_JC,2));
		for ii = 1:size(selected_twins_JC,2)
			act = [];
			iwin = selected_twins_JC(ii);
			twin = timewindows_JC{iwin}{1};
			eval(['act = [act; tuning.JC.',pairname_JC,'.neuract.bytrial.',twin,'];'])
			fr_JC(ii) = mean(act(:,end));
		end
        %SO trials
        fr_SO = nan(1,size(selected_twins_SO,2));
		for ii = 1:size(selected_twins_SO,2)
			act = [];
			iwin = selected_twins_SO(ii);
			twin = timewindows_SO{iwin}{1};
			eval(['act = [act; tuning.SO.',pairname_SO,'.neuract.bytrial.',twin,'];'])
			fr_SO(ii) = mean(act(:,end));
        end
        fr = [fr_JC;fr_SO];
		infocell(ncells).fr = fr;
		
	end %for icell
end %for isession

%output
specs.ncells = ncells;
twins = {};
ii = 0;
for iwin = selected_twins_JC
	ii = ii+1;
	twins(1,ii) = timewindows_JC{iwin}(1);
end
ii = 0;
for iwin = selected_twins_SO
	ii = ii+1;
	twins(2,ii) = timewindows_SO{iwin}(1);
end
specs.twins = twins;
specs.impose_samesign = impose_samesign;

JCSOcellist.infocell = infocell;
JCSOcellist.specs = specs;

toc

%%
%save file
eval(['save ',filesavename,' JCcellist SOcellist JCSOcellist'])
%


