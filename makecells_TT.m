function makecells_TT(varargin)
% makecells_TT
%
% reads the monkeylogic and files and the spike2 files and saves for each cell
% only one file ***data.mat file, which contains all the cell data for the session
% in the variables sessionParams, psyphydata, trialRecord and celldata

%
% author: camillo, 01/2004.
% revised 11/2008, 12/2010,
% revised by SB: 12/2016
% revised by WS: 10/2018 for TT
% fixed some problems with Juan TT recording by WS: 08/2019

if ~isempty(varargin)
	sessions = varargin;
    
else
	clear all %#ok<CLALL>
	sessions = {
		'J190731a'
		};
end


for isession = 1:size(sessions,1)
	session = sessions{isession};
	readsession_TT
	disp(' ')
	disp(['... SESSION ',filename])
	
	%load ML files, bhv
	disp('... reading behavioral files')
	filename = [dirroot, fileroot_ML, '.bhv2'];
    [bhv,~,TrialRecord,~] = mlread(filename);
    %bhv = bhv_read(filename);
	%bhv = bhv_read_rig2(filename);
	%bhv = bhv_read_rig2bis(filename);
    
    % fix the problem: first trial is missed in TrialRecord compared with
    % bhv %% WS
    for i = 1:length([TrialRecord.User.trialRecord.trialNumber])
        TrialRecord.User.trialRecord(i).trialNumber = TrialRecord.User.trialRecord(i).trialNumber+1;
    end
    
	trials_bhv = [bhv.Trial]';
	ntrials_bhv = size(trials_bhv, 1);
	if ~isequal(trials_bhv, (1:ntrials_bhv)')
		disp('some problem with trial numbers in .bhv file')
		dummy
    end
    
	abstime= [[bhv.AbsoluteTrialStartTime]' [bhv.TrialError]'];
	
	%
	%unpack codeTimes, codeNumbers
	behav = [];
    codeTimes = [];
    codeNumbers = [];
	for trial = 1:ntrials_bhv
        codeTimes{trial} = bhv(trial).BehavioralCodes.CodeTimes;
        codeNumbers{trial} = bhv(trial).BehavioralCodes.CodeNumbers;
		behav = [behav; codeTimes{trial}, trial*ones(size(codeTimes{trial})), codeNumbers{trial}]; %#ok<AGROW>
	end
% 	%remove redundant codes 9 and 18
% 	for redundant_code = [9, 18]
% 		aux = behav(:,3)==redundant_code;
% 		ind = aux .* [aux(2:end);aux(1)] .* [aux(3:end);aux(1:2)];
% 		jnd = [0; ind(1:end-1)] + [0; 0; ind(1:end-2)];	%1 for redundant codes
% 		behav = behav(~jnd,:);
% 		if ~isequal(sum(aux),3*sum(ind))
% 			disp(['some problem with codes ', num2str(redundant_code)])
% 		end
% 	end

	%make time relative to trial start
	ind = find(behav(:,3)==9);
	for i = 1:length(ind)
		ii = ind(i);
		try jj = ind(i+1); catch, jj = size(behav,1) + 1; end %#ok<*CTCH>
		t0 = behav(ii,1);
		behav(ii:jj-1, 1) = behav(ii:jj-1, 1) - t0;
	end
	behav_bhv = behav;
	
    cycleRate = [];
	for trial = 1:ntrials_bhv
		cycleRate(trial) = bhv(trial).CycleRate(2); %bug in monkeylogic, minCycleRate is actually average cycle rate
	end
	trialError = [bhv.TrialError];
	
	%load ML files, mat
% 	filename = [dirroot, fileroot_ML, '.mat'];
% 	eval(['load ',filename])
	sessionParams = TrialRecord.User.sessionParams;
	trialRecord_mat = TrialRecord.User.trialRecord;
	trials_mat = [trialRecord_mat.trialNumber]';
	if size(trials_mat,1)~=ntrials_bhv &&  size(trials_mat,1)~=ntrials_bhv+1 &&  size(trials_mat,1)~=ntrials_bhv-1
		disp('mismatch with number of trials between .bhv and .mat files')
	end
	
	%load CED files
	disp('... reading spike files')
	%initialize
	behav_CED = [];
	ncells = size(cells_td, 1); %#ok<NODEF>
	celldata_all = cell(1,ncells);
	for icell = 1:ncells
		celldata_all{icell} = [];
	end
	
	if isfield(parsession,'no_CED')
    else
        for ifile = filenums_CED
            if	ifile<10,	ifilestr = ['0', num2str(ifile)];
            else			ifilestr = num2str(ifile);
            end
				
            %	ifilestr = num2str(ifile);
            filename = [dirroot, fileroot_CED, ifilestr, '.mat'];
            eval(['load ',filename])
				
            %behav
            try
                varname = [fileroot_CED, ifilestr, '_Ch32'];
                eval(['markerchannel = ', varname, ';'])
            % 		catch
            % 			varname = [fileroot_CED, ifilestr, 'fix_Ch32'];
            % 			eval(['markerchannel = ', varname, ';'])
            catch
                varname = [fileroot_CED, ifilestr, '_32_bit__Ch32'];
                eval(['markerchannel = ', varname, ';'])
            end
            behav = [markerchannel.times double(markerchannel.codes(:,1))];
            %

%             %LFP
%             try
%                 for k=1:4
%                 varname = [fileroot_CED, ifilestr, '_Ch4' num2str(k)];
%                 eval(['lfpchannel = ', varname, ';'])
%                 LFP{k}= lfpchannel.values;
%                 end
%             catch
%                 LFP{k}= [];
%             end
%             filename = [dirroot, session, '_LFP'];
%             eval(['save ', filename , ' LFP behav'])
		
            %data = {}; mtrace={}; minth={}; minfo={}; mval={};
            for icell = 1:ncells
                cellId = cells_td(icell, :);
                try
                    if ~uprobe
                        varname = [fileroot_CED, ifilestr, '_Ch2', num2str(cellId(1))];
                        eval(['elchannel = ', varname, ';'])
                    elseif uprobe
                        varname = [fileroot_CED, ifilestr, '_Ch', num2str(10+cellId(1))];
                        eval(['elchannel = ', varname, ';'])
                    end
                    % 			catch
                    % 				varname = [fileroot_CED, ifilestr, 'fix_Ch2', num2str(cellId(1))];
                    % 				eval(['elchannel = ', varname, ';'])
                catch
                    varname = [fileroot_CED, ifilestr, '_32_bit__Ch', num2str(30+cellId(1))];
                    eval(['elchannel = ', varname, ';'])
                end
			
                %keyboard
                eldata = [elchannel.times double(elchannel.codes(:,1))];
                ind = find(eldata(:,2) == cellId(2));
                data{icell} = [eldata(ind,1), cellId(1)*ones(size(ind,1),1), eldata(ind,2)]; %#ok<AGROW>
			
                		
            end % for icell
        
			
			%keyboard
			%
			[behav2, datacell] = ced_psyphydecode(behav, data);
			behav_CED = [behav_CED; behav2]; %#ok<AGROW>
% 			behav_CED = behav2; 
			for icell = 1:ncells
				celldata_all{icell} = [celldata_all{icell}; datacell{icell}];
			end
        end
        
        % fix problematic sessions of Juan_TT
        if isequal(session,'J190731a')
            ind_fix = find(behav_CED(:,2)<26); 
            behav_CED(ind_fix,:) = []; % trial 1-26 is missed in bhv
            behav_CED(:,2) = behav_CED(:,2)-25; % trial 1-26 is missed in bhv
        end
        if isequal(session,'J190801c')
           ind_fix = find(behav_CED(:,2)>466); 
           behav_CED(ind_fix,2) = behav_CED(ind_fix,2)+1; % trial 467 is missed in CED
        end
        trials_CED = unique(behav_CED(:,2));       
	end
		
	%
	%find good trials
    %remove trials with slow cycle rate
	min_cycleRate = 1000;	%min acceptable cycle rate
	ind_fast = find(cycleRate >= min_cycleRate);
    if length(ind_fast)/ntrials_bhv < .9
		disp('note: more than 10% trials removed because of slow cycle rate')
    end
	trials_bhv = trials_bhv(ind_fast);
    trialError = trialError(ind_fast);
	ind = ismember(behav_bhv(:,2), trials_bhv);
    behav_bhv = behav_bhv(ind,:);
	%	
    %remove errors
	remove_errors = 1;
    if remove_errors
		ind_hits = ~trialError;
        trials_bhv = trials_bhv(ind_hits);
		ind = ismember(behav_bhv(:,2),trials_bhv);
		behav_bhv = behav_bhv(ind,:);
    end
	
	%
    %intersect bhv, mat and CED data	
    if isfield(parsession,'no_CED')
	    trials =trials_bhv;
    else
        %
        %intersect bhv, mat and CED data
        trials = intersect(trials_CED, intersect(trials_mat, trials_bhv));
        ind = ismember(behav_bhv(:,2),trials);
        behav_bhv = behav_bhv(ind,:);
        ind = ismember(behav_CED(:,2),trials);
        behav_CED = behav_CED(ind,:);
                  
        if ~isequal(size(behav_bhv,1), size(behav_CED,1))
            badtrials = [];
            % 		for itrial = 1:length(trials)
            for itrial = trials'
                ind = behav_bhv(:,2)==itrial;
                jnd = behav_CED(:,2)==itrial;
                if ~isequal(sum(ind),sum(jnd))
                  	badtrials = [badtrials; itrial]; %#ok<AGROW>
                end
            end
            trials = setdiff(trials, badtrials);
            ind = ismember(behav_bhv(:,2),trials);
            behav_bhv = behav_bhv(ind,:);
            ind = ismember(behav_CED(:,2),trials);
            behav_CED = behav_CED(ind,:);
            disp(['removed ',num2str(length(badtrials)),' trials because of mismatch between CED and ML files'])
        end
        %
        %check that time correspond between bhv and CED
        max_tdiff = 5;	%(ms) max acceptable time mismatch between bhv and CED
        tdiff = behav_bhv(:,1) - behav_CED(:,1);
        ind = abs(tdiff)>max_tdiff;
        notr = unique(behav_bhv(ind,2));
        if ~isempty(notr)
            trials = setdiff(trials, notr);
            ind = ismember(behav_bhv(:,2),trials);
            behav_bhv = behav_bhv(ind,:);
            ind = ismember(behav_CED(:,2),trials);
            behav_CED = behav_CED(ind,:); %#ok<NASGU>
        end
    end
    psyphydata = behav_CED;			%should it be behav_CED?? % SB bhV>CED
    
    
				
	%fix session?
    try eval(['fixsession_', session]), catch, end
		
	%
    %compute trialRecord, goodTrials and celldata and save datafiles
	disp('... saving cell data files')
    trials_aux = trials;	%makes sure that fixcell does not affect subsequent cells
    
    
    
	psyphydata_aux = psyphydata;
    for icell = 1:ncells
		cellId = cells_td(icell,:);
        cellname = [num2str(cellId(1)), num2str(cellId(2))];
		%if necessary remove initial trials (see parsession)
		if isfield(parsession,'rem_firsttrials')
            trials = trials(parsession.rem_firsttrials:end);
		end
		%fix cell?
		try eval(['fixcell_', session, cellname]);  catch, end
        
                                
		ind = ismember(trials_mat, trials);
		trialRecord = trialRecord_mat(ind); 
			
% 		% Create trial error % SB
% 		trialRecordError = trialRecord_mat(trialError>0);
% 		trialerrors=trialError(trialError>0);
% 		for k=1:size(trialRecordError,2)
% 		[trialRecordError(k).errortype]=deal(trialerrors(k));
% 		end
%		
			
        % seperate JC and SO
        goodTrials_JC = TrialRecord.User.goodTrials_JC; %keyboard
        if size(goodTrials_JC,2) == 6
            goodTrials_JC(:,4) = [];
        end
        goodTrials_JC(:,1) = goodTrials_JC(:,1)+1; % miss the first trial in TrialRecord compared with in bhv data %WS
        % also fix in line 40
        goodTrials_SO = TrialRecord.User.goodTrials_SO; %keyboard
        if size(goodTrials_SO,2) == 6
            goodTrials_SO(:,4) = [];
        end
        goodTrials_SO(:,1) = goodTrials_SO(:,1)+1; %% miss the first trial in TrialRecord compared with in bhv data %WS
        % also fix in line 40
    	jnd_JC = find(ismember(goodTrials_JC(:,1), trials));
        jnd_SO = find(ismember(goodTrials_SO(:,1), trials));
        % if ~isequal(size(jnd_JC,1)+size(jnd_SO,1), size(trials,1))
        if (size(jnd_JC,1)+size(jnd_SO,1)-size(trials,1))<-1 % %% miss the first trial in TrialRecord compared with in bhv data %WS
			disp('some problem with goodTrials_JC (too few!)')
			dummy
        end
		goodTrials_JC = goodTrials_JC(jnd_JC,:);
		goodTrials_SO = goodTrials_SO(jnd_SO,:);
        %
		if isfield(parsession,'no_CED')
			celldata=[]; celldataerror =[];
        else
            ind = ismember(celldata_all{icell}(:,2), trials);
            celldata = celldata_all{icell}(ind,:);
            celldataerror = celldata_all{icell}(~ind,:);
            %
            ind = ismember(psyphydata(:,2), trials);
            psyphydata = psyphydata(ind,:); %#ok<NASGU>
            %
        end
    	filename = [dirroot, session, cellname, '_data'];
			
			
% 		eyesignal=[];
% 		%%%% ADD EYES TRACES
% 		for i=1:numel(trials)			
% 			ind=find(codeNumbers{trials(i)}==39);
% 			eyebegin=codeTimes{trials(i)}(ind);
% 			eyesign=bhv.AnalogData{trials(i)}.EyeSignal((eyebegin-1000):(eyebegin+500),:);
% 			eyesignal(i).X=eyesign(:,1);
% 			eyesignal(i).Y=eyesign(:,2);
% 		end
			
		% trialserror=find(trialError>0);
%       for i=1:numel(trialserror)
%       trialRecordError(i).eyes=single(bhv.AnalogData{trialserror(i)}.EyeSignal);
%       end
			
%       keyboard
        if exist('trace')
            trace=[];
        end
        
        % remove extra juice(50) in psyphydata
        for nn = size(psyphydata,1):-1:2
            if psyphydata(nn,3)==psyphydata(nn-1,3)
                psyphydata(nn,:)=[];
            end
        end
                
		% eval(['save ', filename , ' eyesignal psyphydata goodTrials_JC goodTrials_SO trialRecord celldata trialRecordError celldataerror trace'])
        eval(['save ', filename , ' psyphydata goodTrials_JC goodTrials_SO trialRecord celldata celldataerror trace'])
		%
        % get table01_all, pairnames and save file *_bhvParams
        % 		[goodtable goodlabels juices] = get_goodtable(sessionParams);
		%	keyboard
		[~, goodlabels, ~] = get_goodtable(sessionParams);
            
        % seperate JC and SO
        [pairlist_JC, table01_all_JC] = make_choicepattern_TT(goodTrials_JC,'JC'); %#ok<ASGLU>
		npairs_JC = size(pairlist_JC, 1);
		pairs_JC = cell(npairs_JC,1);
		for ii = 1:npairs_JC
            pairs_JC{ii} = [goodlabels{pairlist_JC(ii,1)}(1),goodlabels{pairlist_JC(ii,2)}(1)];
        end
        [pairlist_SO, table01_all_SO] = make_choicepattern_TT(goodTrials_SO,'SO'); %#ok<ASGLU>
        npairs_SO = size(pairlist_SO, 1);
		pairs_SO = cell(npairs_SO,1);
		for ii = 1:npairs_SO
            pairs_SO{ii} = [goodlabels{pairlist_SO(ii,1)}(1),goodlabels{pairlist_SO(ii,2)}(1)];
		end
		filename = [dirroot, session, cellname, '_bhvParams'];
		eval(['save ', filename , ' sessionParams pairs_JC pairlist_JC table01_all_JC pairs_SO pairlist_SO table01_all_SO'])	
        %
		trials = trials_aux;	%resets trials and psyphydata
		psyphydata = psyphydata_aux;
			
    end %for icell
	
    disp(['... saved ',num2str(ncells),' cell data files for session: ',session]);
	
end %for isession
		
	
	
