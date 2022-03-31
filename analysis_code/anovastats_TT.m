function [anovastats] = anovastats_TT(cellname, atleast_nntrials) %#ok<STOUT>
%
% this function makes the anovastats

% author: camillo, september 2004.
% revisions january 2005, 07, 08, 09, 10
% 			april 2013		adapted for RM
%           December 2016	adapted for SO
%           October 2018    adapted for TT % WS

%  keyboard

if 0
	clear all %#ok<*UNRCH>
	cellname = 'V101209a11';
	atleast_nntrials = 2;
end

disp(['   ... computing anovastats for cell ',cellname])
%
JCSO = {'SO','JC'};

for iJCorSO = 1:length(JCSO)
    JCorSO = JCSO{iJCorSO};
    
    [timewindows, ~] = get_timewindows_TT(JCorSO);
    ntwins = length(timewindows);

	%load sigmoidfit, tuning
	session = cellname(1:8); readsession_TT 
	filename = [dirroot,cellname,'_tuning']; eval(['load ',filename])
    if isequal(JCorSO,'JC')
        pairnames = {'AB'};
    elseif isequal(JCorSO,'SO')
        pairnames = {'AB','BA','ABA'}; %fieldnames(tuning);
    end
	%if length(pairname)>1,dummy, end
	

	%analyses to do, nfactors
    if isequal(JCorSO,'JC')
        ana_todo = {
            'byoffertypeNside',		3;
            'bytrialtype',			1;
            };
    elseif isequal(JCorSO,'SO')
        ana_todo = {
            %'byOffertypeOrderSide',		6;
            %'byChoosenJuiceOrderSide',		7;
            'byOffertypeOrderSideChoosen',		4;
            'byOfOrChinteractions', 4
            'bytrialtype', 1;
            };
    end
        
	
    for ipair= 1:numel(pairnames)	
    %initialize
        for iana = 1:size(ana_todo,1)
            currana = ana_todo{iana,1};
            nfactors = ana_todo{iana,2}; %#ok<NASGU>
            eval(['anovastats.',JCorSO,'.',pairnames{ipair},'.',currana,'.tab = [];'])
            eval(['anovastats.',JCorSO,'.',pairnames{ipair},'.',currana,'.pval = nan(ntwins,nfactors);'])
        end


        for iwin = 1:ntwins
            twin = timewindows{iwin}{1};
	
            %get neuract and remove forced choices redundancy
            % colomn in neuract.bytrial:
            % JC: trial# #A #B posA chosenpos order(only=1) FR
            % SO: trial# #A #B posA chosenID orderA FR
            if ~isequal(cellname(8),'x')
                eval(['neuract = tuning.',JCorSO,'.',pairnames{ipair},'.neuract.bytrial.',twin,';'])
                [~,ii] = unique(neuract(:,1));
                neuract = neuract(sort(ii),:);
%               neuract(~neuract(:,4),2) = 0;			%convention for forced choices
%               neuract(~neuract(:,5),3) = 0;
            else
                eval(['neuract = tuning_aux.XX.neuract.bytrial.',twin,';'])
                [~,ii] = unique(neuract(:,1:2),'rows');
                neuract = neuract(sort(ii),:);
                neuract(~neuract(:,5),3) = 0;			%convention for forced choices
                neuract(~neuract(:,6),4) = 0;
		%
                eval(['neuract_norm = tuning_aux.XX.neuract_norm.bytrial.',twin,';'])
                [~,ii] = unique(neuract_norm(:,1:2),'rows');
                neuract_norm = neuract_norm(sort(ii),:);
                neuract_norm(~neuract_norm(:,5),3) = 0;	%convention for forced choices
                neuract_norm(~neuract_norm(:,6),4) = 0;
            end
	
            for iana = 1:size(ana_todo,1)
                currana = ana_todo{iana,1};
                nfactors = ana_todo{iana,2};
		
                if isequal(JCorSO,'JC')
                    if		isequal(currana,'byoffertypeNside'),	unitypes_all = neuract(:,2:5);
                    elseif	isequal(currana,'bytrialtype'),			unitypes_all = [neuract(:,2:3), prod(neuract(:,4:5),2)];
                    end
                elseif isequal(JCorSO,'SO')
                    %if		isequal(currana,'byOffertypeOrderSide'),		unitypes_all = [neuract(:,2:3), neuract(:,4) neuract(:,6)];
                    %elseif	isequal(currana,'byChoosenJuiceOrderSide'),		unitypes_all = [neuract(:,2:3) neuract(:,4) prod(neuract(:,5:6),2) prod(neuract(:,4:5),2) ];		
                    if		isequal(currana,'byOffertypeOrderSideChoosen'),		unitypes_all = [neuract(:,2:6)];
                    elseif	isequal(currana,'byOfOrChinteractions'),		unitypes_all = [neuract(:,2:3), neuract(:,5:6) ];
                    elseif	isequal(currana,'bytrialtype'),			unitypes_all = [neuract(:,2:3), neuract(:,5:6) ];

                    %elseif	isequal(currana,'byoffertypeXNside'),	unitypes_all = neuract(:,[1,3:6]);
                    %elseif	isequal(currana,'bytrialtypeX'),		unitypes_all = [neuract(:,[1,3:4]), prod(neuract(:,5:6),2)];
                    %elseif	isequal(currana,'bytrialtypeX_norm'),	unitypes_all = [neuract_norm(:,[1,3:4]), prod(neuract_norm(:,5:6),2)];
                    %elseif	isequal(currana,'bysession'),			unitypes_all = neuract(:,1);
                    end
                end
    
                
		
                %apply atleast_nntrials criterion
                [~,~,typenums_all] = unique(unitypes_all,'rows');
                ntypes = size(unique(typenums_all),1);
                mult = zeros(ntypes,1);
                for itype = 1:ntypes, mult(itype) = sum(typenums_all==itype); end
                typenums_good = find(mult>atleast_nntrials);
                [ind, loc] = ismember(typenums_all, typenums_good);
		
                if sum(ind)==0
                    disp('NOT ENOUGHT TRIALS!!!!!!!')
                    % keyboard
                end
		
                if	~isequal(currana,'bytrialtypeX_norm')
                    neuract_aux = neuract(ind,:);
                else
                    neuract_aux = neuract_norm(ind,:);
                end
                act = neuract_aux(:,end);

                %n-way anova
                if sum(act>0)>1
                    if ~isequal(cellname(8),'x')
                        [~, ~, offtypenum] = unique(neuract_aux(:,2:3),'rows');
                        if isequal(JCorSO,'JC')
                            pos = neuract_aux(:,4);
                            mov = neuract_aux(:,5);
                            trialtypenum = typenums_all(ind,:);
                        elseif isequal(JCorSO,'SO')
                            pos = neuract_aux(:,4);
                            chosAB= neuract_aux(:,5);
                            ord= neuract_aux(:,6);
                            chospos = neuract_aux(:,5).*neuract_aux(:,4);
                            chosord= neuract_aux(:,5).*neuract_aux(:,6);
                            posord= neuract_aux(:,4).*neuract_aux(:,6);
                            trialtypenum = typenums_all(ind,:);
                        end
                    else
                        [~, ~, offtypeXnum] = unique(neuract_aux(:,[1,3:4]),'rows');
                        pos = neuract_aux(:,4);
                        mov = neuract_aux(:,5).*neuract_aux(:,4);
                        trialtypeXnum = typenums_all(ind,:);
                        sessnum = unitypes_all(ind,1);
                    end
			
                    if isequal(JCorSO,'JC')
                        if		isequal(currana,'byoffertypeNside')
                            [pval, tab, stats, terms] = anovan(act, {offtypenum, pos, mov}, 'linear', 3, {'offertype','pos','mov'}, 'off'); %#ok<*ASGLU>			
			
                        elseif	isequal(currana,'bytrialtype')
                          %keyboard
                            [pval, tab, stats, terms] = anovan(act, {trialtypenum},'linear', 3, {'trialtype'}, 'off');
                        end

                    elseif isequal(JCorSO,'SO')
%                       if		isequal(currana,'byOffertypeOrderSide') && ipair==3
%                       [pval, tab, stats, terms] = anovan(act, {offtypenum,ord,pos}, ...
%                       'interaction',3, {'Offertype', 'Order', 'Side'}, 'off'); %#ok<*ASGLU>
% 	
%                       elseif	isequal(currana,'byChoosenJuiceOrderSide') && ipair==3
%                           [pval, tab, stats, terms] = anovan(act, {chosAB,chosord,chospos}, ...
%                           'full',3, {'ChJuice','ChOrder', 'ChSide'}, 'off'); %#ok<*ASGLU>
		
                        if		isequal(currana,'byOffertypeOrderSideChoosen') %&& ipair==3                        
                            [pval, tab, stats, terms] = anovan(act, {offtypenum,chosord,chosAB,chospos}, ...
                            'linear',3,{'Offertype', 'ChOrder','ChJuice','ChSide'}, 'off'); %#ok<*ASGLU>
	
                        elseif	isequal(currana,'byOfOrChinteractions') %&& ipair==3			
                            [pval, tab, stats, terms] = anovan(act, {offtypenum,ord,chosAB}, ...
                            [1 1 0; 1 0 1; 0 1 1; 1 1 1], 2 ,{'Offertype','Order','ChJuice'}, 'off'); %#ok<*ASGLU>
			
                        elseif	isequal(currana,'bytrialtype')			
                            [pval, tab, stats, terms] = anovan(act, {trialtypenum}, ...
                            'linear', 3, {'trialtype'}, 'off');
				
                        else
                            pval = NaN*ones(nfactors,1); %#ok<NASGU>
                            tab = NaN*ones(nfactors,1); %#ok<NASGU>
                        end	
                    end
			
                else
                    pval = NaN*ones(nfactors,1); %#ok<NASGU>
                    tab = NaN*ones(nfactors,1); %#ok<NASGU>
                end

                %output
                try		
                    eval(['anovastats.',JCorSO,'.',pairnames{ipair},'.',currana,'.tab.',timewindows{iwin}{1},' = tab;'])
                    eval(['anovastats.',JCorSO,'.',pairnames{ipair},'.',currana,'.pval(iwin,:) = pval;'])
                catch
                    keyboard
                end
		
            end %for iana
        end %for twin
    end %for ipair
end %for iJCorSO