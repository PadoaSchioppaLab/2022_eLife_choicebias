function [] = tuning_TT(cellname)
%
% makes and saves the variable and tuning (including activity trial by trial).
% several tuning curves are computed: 
% byoffertype:	distinguishes offertypes
% bytrialtype:	distinguishes A/B choices -- most of the analysis based on this.
% by trial:		trial by trial.

%
% author: camillo, january 2004
% revision: september 2004
% 			september 2005
%			december 2008: adapted for RC
%			december 2009: adapted for DT
%			march 2010: revised for radial analysis
%			december 2010: 'downgraded' to JC
%           december 2016 adapted to SO
%           october 2018 adapted to TT
%

% keyboard

if 0
	clear all %#ok<UNRCH>
	cellname = 'G180910a11';
end

disp(['   ... computing tuning of cell ',cellname])

session = cellname(1:8); readsession_TT 
% [timewindows colrs] = get_timewindows(cellname);
%
% 
JCSO = {'JC','SO'};

for iJCorSO = 1:length(JCSO)
    JCorSO = JCSO{iJCorSO};
    
    [timewindows, ~, timewindows_off, timewindows_cho, ~, ~, timewindows_bas] = get_timewindows_TT(JCorSO);
    ntwins = length(timewindows);
    ntwins2 = length(timewindows_off);
    ntwins3 = length(timewindows_cho);
%     ntwins4 = length(timewindows_tgt);
%     ntwins5 = length(timewindows_cue);
    ntwins6 = length(timewindows_bas);
    
    %load cell data, bhvParams
    filename = [dirroot, cellname, '_data']; eval(['load ',filename])
    filename = [dirroot, cellname, '_bhvParams']; eval(['load ',filename])
    ngoods = size(sessionParams.goods,2);
    
    if isequal(JCorSO,'JC')
        npairs = 1; 
        pairs={'AB'};
    elseif isequal(JCorSO,'SO')
        npairs = 3; 
        pairs={'BA' 'AB' 'ABA'};
        % pairs={'AB' 'BA' 'ABA'};
    end
    
    %in psyphydata convert possible_outcomes to flag_outcome
    ind = ismember(psyphydata(:,3), flags.poss_outcomes);
    psyphydata(ind,3) = flags.outcome;

    for ipair = 1:npairs
        if ipair ==1
            if isequal(JCorSO,'JC')
                order = 0; 
            elseif isequal(JCorSO,'SO')
                order = 1;  % Select BA by: hitdata(:,6)~=order
            end
        elseif ipair==2 
             order=-1; % Select AB by: hitdata(:,6)~=order
        elseif ipair==3; order=0;  end % SELECT both by: hitdata(:,6)~=order
		
        if order~=0 % if ipair~=3
            eval(['table01 = abs(table01_all_',JCorSO,'.bydir{-(ipair-3)});']) % normal order is first AB then BA, HERE WE DO OPPOSITE
            eval(['table02= table01_all_',JCorSO,'.bydir{-(ipair-3)};'])
        else
            eval(['table01 = table01_all_',JCorSO,'.pooldir{1};'])
            eval(['table02= vertcat((table01_all_',JCorSO,'.bydir{1}),(table01_all_',JCorSO,'.bydir{2}));'])
        end
	
        eval(['tuning.',JCorSO,'.',pairs{ipair},'.table01 = table01;']) %{ipair};
        eval(['tuning.',JCorSO,'.',pairs{ipair},'.table02 = table02;']) %{ipair};
        eval(['hitdata = get_hitdata_TT(goodTrials_',JCorSO,', pairlist_',JCorSO,'(1,:), ngoods, JCorSO);']);
        if isequal(JCorSO,'JC')
            hitdata(:,6)=1; 
        end
        
        % regular tuning analysis, different time windows
        for iwin = 1 : ntwins
            twin = timewindows{iwin}{1};
            activity = []; 
            tuncurve_offer = [];
            tuncurve_trial = [];
		
            for iofftype = 1 : size(table01,1)
                %tuncurve by offer type
                ind = find(hitdata(:,2)==table01(iofftype,1) & hitdata(:,3)==table01(iofftype,2) & hitdata(:,6)~=order); 
                if ~isempty(ind)
                    hits = hitdata(ind,1);
                    act = actmake_TT(celldata, psyphydata, hits, timewindows{iwin});
                    tun = [table01(iofftype,[1,2]), mean(act), std(act), length(hits)];
                    tuncurve_offer = [tuncurve_offer; tun];
                else
                    disp('empty ind')
                end
			
                %%%% 
                if order~=0 % ipair~=3 % for BA and AB in SO
                    %tuncurve by trial type
                    for itaste = [-1,1] % -1 for B % WS
                        ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                            hitdata(:,3)==table01(iofftype,2) & hitdata(:,5)==itaste & hitdata(:,6)~=order );
                        if ~isempty(ind)
                            hits = hitdata(ind,1);
                            act = actmake_TT(celldata, psyphydata, hits, timewindows{iwin});
                            tun = [hitdata(ind(1),2:3), itaste, order*-1, mean(act), std(act), length(hits)];
                            tuncurve_trial = [tuncurve_trial; tun];
                            %
                            %activity, trial by trial
                            act = actmake_TT(celldata, psyphydata, hits, timewindows{iwin});
                            typeact = [hitdata(ind,:),act];
                            activity = [activity; typeact];			
                        end
                    end %for itaste
			
                else % for JC and ABA(both) in SO
                    %tuncurve by trial type
                    for itaste = [-1,1]
                        if isequal(JCorSO, 'JC')
                            ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                                hitdata(:,3)==table01(iofftype,2) & ...
                                hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)~=order );
                            if ~isempty(ind)
                                hits = hitdata(ind,1);
                                act = actmake_TT(celldata, psyphydata, hits, timewindows{iwin});
                                tun = [hitdata(ind(1),2:3), itaste, order, mean(act), std(act), length(hits)];
                                tuncurve_trial = [tuncurve_trial; tun];
                                %
                                %activity, trial by trial
                                act = actmake_TT(celldata, psyphydata, hits, timewindows{iwin});
                                typeact = [hitdata(ind,:),act];
                                activity = [activity; typeact];
                            end				
                        elseif isequal(JCorSO, 'SO')
                            for ord= [-1,1]
                                ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                                hitdata(:,3)==table01(iofftype,2) & ...
                                hitdata(:,5)==itaste & hitdata(:,6)==ord );%%???
                                %hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
                                if ~isempty(ind)
                                    hits = hitdata(ind,1);
                                    act = actmake_TT(celldata, psyphydata, hits, timewindows{iwin});
                                    tun = [hitdata(ind(1),2:3), itaste, ord, mean(act), std(act), length(hits)];
                                    tuncurve_trial = [tuncurve_trial; tun];
                                    %
                                    %activity, trial by trial
                                    act = actmake_TT(celldata, psyphydata, hits, timewindows{iwin});
                                    typeact = [hitdata(ind,:),act];
                                    activity = [activity; typeact];
                                end				
                            end %for ord
                        end % if isequal(JCorSO, 'JC')
                    end %for itaste			
                end			
            end %for iofftype
		
		%output
		try
            eval(['tuning.',JCorSO,'.',pairs{ipair},'.neuract.bytrial.',			twin,' = activity;'])
            eval(['tuning.',JCorSO,'.',pairs{ipair},'.neuract.byoffertype.',		twin,' = tuncurve_offer;'])
            eval(['tuning.',JCorSO,'.',pairs{ipair},'.neuract.bytrialtype.',		twin,' = tuncurve_trial;'])
		catch
			keyboard
		end
	end	%for iwin

	
	if 1
% % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETAILED TIMEWINDOWS _OFF
 % % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
		for iwin = 1 : ntwins2
		twin = timewindows_off{iwin}{1};
		activity = []; 
		tuncurve_offer = [];
		tuncurve_trial = [];
		
		for iofftype = 1 : size(table01,1)
			%tuncurve by offer type
			ind = find(hitdata(:,2)==table01(iofftype,1) & hitdata(:,3)==table01(iofftype,2) & hitdata(:,6)~=order);
			if ~isempty(ind)
				hits = hitdata(ind,1);
				act = actmake_TT(celldata, psyphydata, hits, timewindows_off{iwin});
				tun = [table01(iofftype,[1,2]), mean(act), std(act), length(hits)];
				tuncurve_offer = [tuncurve_offer; tun];
			else
				disp('empty ind')
			end
			
			if order~=0 % ipair~=3 % for BA and AB in SO
			%tuncurve by trial type
			for itaste = [-1,1]
				ind = find(hitdata(:,2)==table01(iofftype,1) & ...
					hitdata(:,3)==table01(iofftype,2) & ...
					hitdata(:,5)==itaste & hitdata(:,6)~=order );
					%hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
				if ~isempty(ind)
					hits = hitdata(ind,1);
					act = actmake_TT(celldata, psyphydata, hits, timewindows_off{iwin});
					tun = [hitdata(ind(1),2:3), itaste, order*-1, mean(act), std(act), length(hits)];
					tuncurve_trial = [tuncurve_trial; tun];
					%
					%activity, trial by trial
					act = actmake_TT(celldata, psyphydata, hits, timewindows_off{iwin});
					typeact = [hitdata(ind,:),act];
					activity = [activity; typeact];			
				end
			end %for itaste	
			
            else
			%tuncurve by trial type
			for itaste = [-1,1]
                if isequal(JCorSO, 'JC')
                    ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                        hitdata(:,3)==table01(iofftype,2) & ...
                        hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)~=order );
                    if ~isempty(ind)
                        hits = hitdata(ind,1);
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_off{iwin});
                        tun = [hitdata(ind(1),2:3), itaste, order, mean(act), std(act), length(hits)];
                        tuncurve_trial = [tuncurve_trial; tun];
                        %
                        %activity, trial by trial
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_off{iwin});
                        typeact = [hitdata(ind,:),act];
                        activity = [activity; typeact];
                    end				
                elseif isequal(JCorSO, 'SO')
                    for ord= [-1,1]
                        ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                            hitdata(:,3)==table01(iofftype,2) & ...
                            hitdata(:,5)==itaste & hitdata(:,6)==ord );
                        %hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
                        if ~isempty(ind)
                            hits = hitdata(ind,1);
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_off{iwin});
                            tun = [hitdata(ind(1),2:3), itaste, ord, mean(act), std(act), length(hits)];
                            tuncurve_trial = [tuncurve_trial; tun];
                            %
                            %activity, trial by trial
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_off{iwin});
                            typeact = [hitdata(ind,:),act];
                            activity = [activity; typeact];
                        end				
                    end %for ord
                end %if isequal(JCorSO, 'JC')
			end %for itaste
			
			end
				
		end %for ioff
		
		%output
		try
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrial.',			twin,' = activity;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.byoffertype.',		twin,' = tuncurve_offer;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrialtype.',		twin,' = tuncurve_trial;'])
		catch
			keyboard
		end
		
		
	end	%for iwin2
	
	end
	
	
	if 0
	
% % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETAILED TIMEWINDOWS _CHO
 % % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for iwin = 1 : ntwins3
		twin = timewindows_cho{iwin}{1};
		activity = []; 
		tuncurve_offer = [];
		tuncurve_trial = [];
		
		for iofftype = 1 : size(table01,1)
			%tuncurve by offer type
			ind = find(hitdata(:,2)==table01(iofftype,1) & hitdata(:,3)==table01(iofftype,2) & hitdata(:,6)~=order);
			if ~isempty(ind)
				hits = hitdata(ind,1);
				act = actmake_TT(celldata, psyphydata, hits, timewindows_cho{iwin});
				tun = [table01(iofftype,[1,2]), mean(act), std(act), length(hits)];
				tuncurve_offer = [tuncurve_offer; tun];
			else
				disp('empty ind')
			end
			
			if order~=0 % ipair~=3 % for BA and AB in SO
			%tuncurve by trial type
			for itaste = [-1,1]
				ind = find(hitdata(:,2)==table01(iofftype,1) & ...
					hitdata(:,3)==table01(iofftype,2) & ...
					hitdata(:,5)==itaste & hitdata(:,6)~=order );
					%hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
				if ~isempty(ind)
					hits = hitdata(ind,1);
					act = actmake_TT(celldata, psyphydata, hits, timewindows_cho{iwin});
					tun = [hitdata(ind(1),2:3), itaste, order*-1, mean(act), std(act), length(hits)];
					tuncurve_trial = [tuncurve_trial; tun];
					%
					%activity, trial by trial
					act = actmake_TT(celldata, psyphydata, hits, timewindows_cho{iwin});
					typeact = [hitdata(ind,:),act];
					activity = [activity; typeact];			
				end
			end %for itaste
			
			
            else
			%tuncurve by trial type
			for itaste = [-1,1]
				if isequal(JCorSO, 'JC')
                    ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                        hitdata(:,3)==table01(iofftype,2) & ...
                        hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)~=order );
                    if ~isempty(ind)
                        hits = hitdata(ind,1);
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_cho{iwin});
                        tun = [hitdata(ind(1),2:3), itaste, order, mean(act), std(act), length(hits)];
                        tuncurve_trial = [tuncurve_trial; tun];
                        %
                        %activity, trial by trial
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_cho{iwin});
                        typeact = [hitdata(ind,:),act];
                        activity = [activity; typeact];
                    end				
                elseif isequal(JCorSO, 'SO')
                    for ord= [-1,1]
                        ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                            hitdata(:,3)==table01(iofftype,2) & ...
                            hitdata(:,5)==itaste & hitdata(:,6)==ord );
                        %hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
                        if ~isempty(ind)
                            hits = hitdata(ind,1);
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_cho{iwin});
                            tun = [hitdata(ind(1),2:3), itaste, ord, mean(act), std(act), length(hits)];
                            tuncurve_trial = [tuncurve_trial; tun];
                            %
                            %activity, trial by trial
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_cho{iwin});
                            typeact = [hitdata(ind,:),act];
                            activity = [activity; typeact];
                        end				
                    end %for ord
                end %if isequal(JCorSO, 'JC')
			end %for itaste
			
			end
			
		end %for ioff
		
		%output
		try
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrial.',			twin,' = activity;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.byoffertype.',		twin,' = tuncurve_offer;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrialtype.',		twin,' = tuncurve_trial;'])
		catch
			keyboard
		end
		
	end	%for iwin3
	
	end
	
	
		if 0
	
% % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETAILED TIMEWINDOWS _TGT
 % % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		for iwin = 1 : ntwins4
		twin = timewindows_tgt{iwin}{1};
		activity = []; 
		tuncurve_offer = [];
		tuncurve_trial = [];
		
		for iofftype = 1 : size(table01,1)
			%tuncurve by offer type
			ind = find(hitdata(:,2)==table01(iofftype,1) & hitdata(:,3)==table01(iofftype,2) & hitdata(:,6)~=order);
			if ~isempty(ind)
				hits = hitdata(ind,1);
				act = actmake_SO(celldata, psyphydata, hits, timewindows_tgt{iwin});
				tun = [table01(iofftype,[1,2]), mean(act), std(act), length(hits)];
				tuncurve_offer = [tuncurve_offer; tun];
			else
				disp('empty ind')
			end
			
			if order~=0 % ipair~=3 % for BA and AB in SO
			%tuncurve by trial type
			for itaste = [-1,1]
				ind = find(hitdata(:,2)==table01(iofftype,1) & ...
					hitdata(:,3)==table01(iofftype,2) & ...
					hitdata(:,5)==itaste & hitdata(:,6)~=order );
					%hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
				if ~isempty(ind)
					hits = hitdata(ind,1);
					act = actmake_TT(celldata, psyphydata, hits, timewindows_tgt{iwin});
					tun = [hitdata(ind(1),2:3), itaste, order*-1, mean(act), std(act), length(hits)];
					tuncurve_trial = [tuncurve_trial; tun];
					%
					%activity, trial by trial
					act = actmake_TT(celldata, psyphydata, hits, timewindows_tgt{iwin});
					typeact = [hitdata(ind,:),act];
					activity = [activity; typeact];			
				end
			end %for itaste
			
			
			else
			
			%tuncurve by trial type
			for itaste = [-1,1]
				if isequal(JCorSO, 'JC')
                    ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                        hitdata(:,3)==table01(iofftype,2) & ...
                        hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)~=order );
                    if ~isempty(ind)
                        hits = hitdata(ind,1);
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_tgt{iwin});
                        tun = [hitdata(ind(1),2:3), itaste, order, mean(act), std(act), length(hits)];
                        tuncurve_trial = [tuncurve_trial; tun];
                        %
                        %activity, trial by trial
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_tgt{iwin});
                        typeact = [hitdata(ind,:),act];
                        activity = [activity; typeact];
                    end				
                elseif isequal(JCorSO, 'SO')
                    for ord= [-1,1] 
                        ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                            hitdata(:,3)==table01(iofftype,2) & ...
                            hitdata(:,5)==itaste & hitdata(:,6)==ord );
                        %hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
                        if ~isempty(ind)
                            hits = hitdata(ind,1);
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_tgt{iwin});
                            tun = [hitdata(ind(1),2:3), itaste, ord, mean(act), std(act), length(hits)];
                            tuncurve_trial = [tuncurve_trial; tun];
                            %
                            %activity, trial by trial
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_tgt{iwin});
                            typeact = [hitdata(ind,:),act];
                            activity = [activity; typeact];
                        end				
                    end %for ord
                end %if isequal(JCorSO, 'JC')
			end %for itaste
			
			end
			
		end %for ioff
		
		%output
		try
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrial.',			twin,' = activity;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.byoffertype.',		twin,' = tuncurve_offer;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrialtype.',		twin,' = tuncurve_trial;'])
		catch
			keyboard
		end
		
	end	%for iwin4
	
    end
	
    
    
    if 0
% % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETAILED TIMEWINDOWS _CUE
 % % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
		for iwin = 1 : ntwins5
		twin = timewindows_cue{iwin}{1};
		activity = []; 
		tuncurve_offer = [];
		tuncurve_trial = [];
		
		for iofftype = 1 : size(table01,1)
			%tuncurve by offer type
			ind = find(hitdata(:,2)==table01(iofftype,1) & hitdata(:,3)==table01(iofftype,2) & hitdata(:,6)~=order);
			if ~isempty(ind)
				hits = hitdata(ind,1);
				act = actmake_TT(celldata, psyphydata, hits, timewindows_cue{iwin});
				tun = [table01(iofftype,[1,2]), mean(act), std(act), length(hits)];
				tuncurve_offer = [tuncurve_offer; tun];
			else
				disp('empty ind')
			end
			
			if order~=0 % ipair~=3 % for BA and AB in SO
			%tuncurve by trial type
			for itaste = [-1,1]
				ind = find(hitdata(:,2)==table01(iofftype,1) & ...
					hitdata(:,3)==table01(iofftype,2) & ...
					hitdata(:,5)==itaste & hitdata(:,6)~=order );
					%hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
				if ~isempty(ind)
					hits = hitdata(ind,1);
					act = actmake_TT(celldata, psyphydata, hits, timewindows_cue{iwin});
					tun = [hitdata(ind(1),2:3), itaste, order*-1, mean(act), std(act), length(hits)];
					tuncurve_trial = [tuncurve_trial; tun];
					%
					%activity, trial by trial
					act = actmake_TT(celldata, psyphydata, hits, timewindows_cue{iwin});
					typeact = [hitdata(ind,:),act];
					activity = [activity; typeact];			
				end
			end %for itaste	
			
            else
			%tuncurve by trial type
			for itaste = [-1,1]
                if isequal(JCorSO, 'JC')
                    ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                        hitdata(:,3)==table01(iofftype,2) & ...
                        hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)~=order );
                    if ~isempty(ind)
                        hits = hitdata(ind,1);
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_cue{iwin});
                        tun = [hitdata(ind(1),2:3), itaste, order, mean(act), std(act), length(hits)];
                        tuncurve_trial = [tuncurve_trial; tun];
                        %
                        %activity, trial by trial
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_cue{iwin});
                        typeact = [hitdata(ind,:),act];
                        activity = [activity; typeact];
                    end				
                elseif isequal(JCorSO, 'SO')
                    for ord= [-1,1]
                        ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                            hitdata(:,3)==table01(iofftype,2) & ...
                            hitdata(:,5)==itaste & hitdata(:,6)==ord );
                        %hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
                        if ~isempty(ind)
                            hits = hitdata(ind,1);
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_cue{iwin});
                            tun = [hitdata(ind(1),2:3), itaste, ord, mean(act), std(act), length(hits)];
                            tuncurve_trial = [tuncurve_trial; tun];
                            %
                            %activity, trial by trial
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_cue{iwin});
                            typeact = [hitdata(ind,:),act];
                            activity = [activity; typeact];
                        end				
                    end %for ord
                end %if isequal(JCorSO, 'JC')
			end %for itaste
			
			end
				
		end %for ioff
		
		%output
		try
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrial.',			twin,' = activity;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.byoffertype.',		twin,' = tuncurve_offer;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrialtype.',		twin,' = tuncurve_trial;'])
		catch
			keyboard
		end
		
		
	end	%for iwin5
    end
    
    if 0
% % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DETAILED TIMEWINDOWS _BAS(ELINE)
 % % % % % % % % % % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
		for iwin = 1 : ntwins6
		twin = timewindows_bas{iwin}{1};
		activity = []; 
		tuncurve_offer = [];
		tuncurve_trial = [];
		
		for iofftype = 1 : size(table01,1)
			%tuncurve by offer type
			ind = find(hitdata(:,2)==table01(iofftype,1) & hitdata(:,3)==table01(iofftype,2) & hitdata(:,6)~=order);
			if ~isempty(ind)
				hits = hitdata(ind,1);
				act = actmake_TT(celldata, psyphydata, hits, timewindows_bas{iwin});
				tun = [table01(iofftype,[1,2]), mean(act), std(act), length(hits)];
				tuncurve_offer = [tuncurve_offer; tun];
			else
				disp('empty ind')
			end
			
			if order~=0 % ipair~=3 % for BA and AB in SO
			%tuncurve by trial type
			for itaste = [-1,1]
				ind = find(hitdata(:,2)==table01(iofftype,1) & ...
					hitdata(:,3)==table01(iofftype,2) & ...
					hitdata(:,5)==itaste & hitdata(:,6)~=order );
					%hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
				if ~isempty(ind)
					hits = hitdata(ind,1);
					act = actmake_TT(celldata, psyphydata, hits, timewindows_bas{iwin});
					tun = [hitdata(ind(1),2:3), itaste, order*-1, mean(act), std(act), length(hits)];
					tuncurve_trial = [tuncurve_trial; tun];
					%
					%activity, trial by trial
					act = actmake_TT(celldata, psyphydata, hits, timewindows_bas{iwin});
					typeact = [hitdata(ind,:),act];
					activity = [activity; typeact];			
				end
			end %for itaste	
			
            else
			%tuncurve by trial type
			for itaste = [-1,1]
                if isequal(JCorSO, 'JC')
                    ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                        hitdata(:,3)==table01(iofftype,2) & ...
                        hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)~=order );
                    if ~isempty(ind)
                        hits = hitdata(ind,1);
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_bas{iwin});
                        tun = [hitdata(ind(1),2:3), itaste, order, mean(act), std(act), length(hits)];
                        tuncurve_trial = [tuncurve_trial; tun];
                        %
                        %activity, trial by trial
                        act = actmake_TT(celldata, psyphydata, hits, timewindows_bas{iwin});
                        typeact = [hitdata(ind,:),act];
                        activity = [activity; typeact];
                    end				
                elseif isequal(JCorSO, 'SO')
                    for ord= [-1,1]
                        ind = find(hitdata(:,2)==table01(iofftype,1) & ...
                            hitdata(:,3)==table01(iofftype,2) & ...
                            hitdata(:,5)==itaste & hitdata(:,6)==ord );
                        %hitdata(:,4).*hitdata(:,5)==itaste & hitdata(:,6)==order );
                        if ~isempty(ind)
                            hits = hitdata(ind,1);
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_bas{iwin});
                            tun = [hitdata(ind(1),2:3), itaste, ord, mean(act), std(act), length(hits)];
                            tuncurve_trial = [tuncurve_trial; tun];
                            %
                            %activity, trial by trial
                            act = actmake_TT(celldata, psyphydata, hits, timewindows_bas{iwin});
                            typeact = [hitdata(ind,:),act];
                            activity = [activity; typeact];
                        end				
                    end %for ord
                end %if isequal(JCorSO, 'JC')
			end %for itaste
			
			end
				
		end %for ioff
		
		%output
		try
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrial.',			twin,' = activity;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.byoffertype.',		twin,' = tuncurve_offer;'])
		eval(['tuning.',JCorSO,'.',pairs{ipair},'.superneuract.bytrialtype.',		twin,' = tuncurve_trial;'])
		catch
			keyboard
		end
		
		
	end	%for iwin6
    end
	
end %for ipair

end %for iJCorSO 
%
%save tuning
filename = [dirroot,cellname,'_tuning_new'];
eval(['save ',filename,' tuning'])

if 0
checksession(session,.8);	%check that the session is complete
end
%tuningplot(cellname);


% 
% %%%%%%  functions  %%%%%%%
% 
% function [act] = actmake(celldata, psyphydata, hits, timewindow)
% %
% % computes the activity of one cell for one particular set of hits 
% % (typically, for one offer type), and for one particular time window.
% %
% act = zeros(length(hits),1);
% 
% for ihit = 1:length(hits)
% 	hit = hits(ihit);
% 	try
% 		ind = psyphydata(:,2)==hit & psyphydata(:,3)==timewindow{2}(1);
% 		timestart = psyphydata(ind,1) + timewindow{3}(1);
% 		ind = psyphydata(:,2)==hit & psyphydata(:,3)==timewindow{2}(2);
% 		timeend = psyphydata(ind,1) + timewindow{3}(2);
% 		
% 		ind = celldata(:,2)==hit;
% 		sts = celldata(ind,:);
% 		
% 		ind = find(sts(:,1)>timestart & sts(:,1)<=timeend);
% 		act(ihit) = 1000*length(ind)/(timeend-timestart);		%activity in sp/sec
% 	catch
% 		disp(['some problem at trial number ',num2str(hit), ' ...see function actmake'])
% 	end
% end


