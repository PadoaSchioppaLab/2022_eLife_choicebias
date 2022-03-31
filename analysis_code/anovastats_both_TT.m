function [anovastats_both] = anovastats_both_TT(cellname, atleast_nntrials) %#ok<STOUT>
%
% this function makes the anovastats

% author: camillo, september 2004.
% revisions january 2005, 07, 08, 09, 10
% 			april 2013		adapted for RM
%           December 2016	adapted for SO
%           October 2018    adapted for TT % WS
%           December 2018   adapted for TT considering both JC and SO  %WS 
%                           only for postoffer and postjuice time windows
%  keyboard

if 0
	clear all %#ok<*UNRCH>
	cellname = 'V101209a11';
	atleast_nntrials = 2;
end

disp(['   ... computing anovastats_both for cell ',cellname])

% 
[timewindows_SO, ~] = get_timewindows_TT('SO');
[timewindows_JC, ~] = get_timewindows_TT('JC');

twins_SO = [2,4,8]; % postoffer1 postoffer2 postjuice
twins_JC = [2,7];   % postoffer postjuice
ntwins = 2; %postoffer and postjuice

ana_todo = { 
             'bytaskNtrialtype',  3;
             'bytrialtype',       1;
             };
    
session = cellname(1:8); readsession_TT;
filename = [dirroot,cellname,'_tuning']; eval(['load ',filename])
tasknames = 'JCSO';    
anovastats_both = [];

% initialize
for iana = 1:size(ana_todo,1)
    currana = ana_todo{iana,1};
    nfactors = ana_todo{iana,2};
    eval(['anovastats_both.',tasknames,'.',currana,'.tab = [];'])
    eval(['anovastats_both.',tasknames,'.',currana,'.pval = nan(ntwins,nfactors);'])
end

for iwin = 1:ntwins
    if iwin == 1    
        twin_JC  = 'postoffer';
        twin_SO1 = 'postoffer1';
        twin_SO2 = 'postoffer2';
        
        %get neuract and remove forced choices redundancy
        % colomn in neuract.bytrial:
        % JC: trial# #A #B posA chosenpos order(only=1) FR
        % SO: trial# #A #B posA chosenID orderA FR
        
        eval(['neuract_JC = tuning.JC.AB.neuract.bytrial.',twin_JC,';'])
        [~,ii] = unique(neuract_JC(:,1));
        neuract_JC = neuract_JC(sort(ii),:);

        eval(['neuract_SO1 = tuning.SO.ABA.neuract.bytrial.',twin_SO1,';'])
        [~,ii] = unique(neuract_SO1(:,1));
        neuract_SO1 = neuract_SO1(sort(ii),:);
        
        eval(['neuract_SO2 = tuning.SO.ABA.neuract.bytrial.',twin_SO2,';'])
        [~,ii] = unique(neuract_SO2(:,1));
        neuract_SO2 = neuract_SO2(sort(ii),:);
        
        neuract_SO = [neuract_SO1(:,1:6),max([neuract_SO1(:,7),neuract_SO2(:,7)]')']; % max value in postoffer1 or postoffer2
        
    elseif iwin == 2
        twin_JC = 'postjuice';
        twin_SO = 'postjuice';
        
        %get neuract and remove forced choices redundancy
        % colomn in neuract.bytrial:
        % JC: trial# #A #B posA chosenpos order(only=1) FR
        % SO: trial# #A #B posA chosenID orderA FR
        
        eval(['neuract_JC = tuning.JC.AB.neuract.bytrial.',twin_JC,';'])
        [~,ii] = unique(neuract_JC(:,1));
        neuract_JC = neuract_JC(sort(ii),:);

        eval(['neuract_SO = tuning.SO.ABA.neuract.bytrial.',twin_SO,';'])
        [~,ii] = unique(neuract_SO(:,1));
        neuract_SO = neuract_SO(sort(ii),:);

    end
    
	
    for iana = 1:size(ana_todo,1)
         currana = ana_todo{iana,1};
         nfactors = ana_todo{iana,2};
		
         if isequal(currana,'bytaskNtrialtype')
             unitypes_all_JC = [neuract_JC(:,2:3), prod(neuract_JC(:,4:5),2), zeros(size(neuract_JC,1),1), -ones(size(neuract_JC,1),1)];  % offerA offerB choice(1 for A) order(0 for JC) tasktype(-1 for JC)
             unitypes_all_SO = [abs(neuract_SO(:,2:3)), neuract_SO(:,5:6), ones(size(neuract_SO,1),1)]; % offerA offerB choice(1 for A) order(1 for AB -1 for BA) tasktype(1 for SO)  
             unitypes_all = [unitypes_all_JC; unitypes_all_SO];       
         elseif  isequal(currana,'bytrialtype')
             unitypes_all_JC = [neuract_JC(:,2:3), prod(neuract_JC(:,4:5),2), zeros(size(neuract_JC,1),1), -ones(size(neuract_JC,1),1)];  % offerA offerB choice(1 for A) order(0 for JC) tasktype(-1 for JC)
             unitypes_all_SO = [abs(neuract_SO(:,2:3)), neuract_SO(:,5:6), ones(size(neuract_SO,1),1)]; % offerA offerB choice(1 for A) order(1 for AB -1 for BA) tasktype(1 for SO)
             unitypes_all = [unitypes_all_JC; unitypes_all_SO];
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
             keyboard
         end
         
         neuract = [neuract_JC; neuract_SO];
         neuract_aux = neuract(ind,:);
         act = neuract_aux(:,end);

         %n-way anova
         if sum(act>0)>1
             trialtypenum = typenums_all(ind,:);
             order01 = unitypes_all(ind,4);
             tasktype = unitypes_all(ind,5);
             [~, ~, offtypenum] = unique(unitypes_all(ind,1:2),'rows');
                        
             
             if isequal(currana,'bytaskNtrialtype')
                 [pval, tab, stats, terms] = anovan(act, {offtypenum, order01, tasktype}, 'linear', 3, {'offertype','order','tasktype'}, 'off'); %#ok<*ASGLU>
                 
             elseif isequal(currana,'bytrialtype')
                 [pval, tab, stats, terms] = anovan(act, {trialtypenum},'linear', 1, {'trialtype'}, 'off');
                 
             end
                   
			
         else
            pval = NaN*ones(nfactors,1); %#ok<NASGU>
            tab = NaN*ones(nfactors,1); %#ok<NASGU>
         end

        %output
        try	
            eval(['anovastats_both.',tasknames,'.',currana,'.tab = [];'])
            eval(['anovastats_both.',tasknames,'.',currana,'.tab.',twin_JC,' = tab;'])
            eval(['anovastats_both.',tasknames,'.',currana,'.pval(iwin,:) = pval;'])
        catch
            keyboard
        end
		
    end %for iana
end %for twin