% % % stumpanalysis_TT
% basic neuron analysis       

close all

tic
monkey = 'both'; % 'Gervinho'; 'Juan'; 'both'
eval(['sessionlist_',monkey]);

%%% brain region: OFC or dLPFC or vLPFC
brain_region = 4; % 0 for all; 4 for OFC; 9 for dLPFC; 10 for vLPFC

%%%

%%%
dosigmoid	 = 1; %OK
dorasters	 = 1; %OK
doprofiles	 = 1; %OK
dotuning	 = 1; %OK
docellstats	 = 1; %OK

%%%
atleast_nntrials = 2; %3 or 2; mostly use 3. Use 2 here because Juan has less trial numbers
mintrialnum = 100;

exnovo		= 0; %OK
savefigs	= 0; %OK
onefile		= 0; %OK
nfiles=0;

fitt='probit';  % Fitting algorithm for sigmoid: 'probit'; OR 'logit';  OR 'normcdf';
waitshow = 0; %OK
limp=0.05; % for ANOVA


%%%%%%
if savefigs && onefile
    filesave = 'C:\Experiments\TwoTasks\Analysis\Analysis_OFC\allcellstats_TT_OFC'; % change based which folder is used 
	filesaveall_save = filesave;
	filesave=[ filesaveall_save '_' num2str(nfiles)];
	eval(['!del ',filesave,'.ps; !del ',filesave,'.pdf']);
end

ncells = 0;
for isession = 1:size(sessions,1)
	
	try
		% 	if (f.bytes*1.1)>=(f2.MaxPossibleArrayBytes)
		f= dir([filesave '.ps']);
		if f.bytes>50000000
			nfiles=nfiles+1;
			filesave=[ filesaveall_save '_' num2str(nfiles)];
		end
	catch
	end
	
	session = sessions{isession};
	disp('	'); disp(['SESSION ',session])
	readsession_TT
	
        
	%if necessary, makecells
% 	try
		if (exnovo), dummy; end %#ok<*UNRCH>
		cellname = [session, num2str(cells_td(1,1) ), num2str(cells_td(1,2))];
		filename = [dirroot, cellname, '_data'];
		eval(['load ',filename])
% 	catch %#ok<*CTCH>
% 		makecells_TT(session);
% 	end
	
    
    % choose cells from different brain regions
    if brain_region == 0
        cells_td = cells_td;
    else
        if ~uprobe
            ind = find(channels_label==brain_region);
            channels_specregion = channels(ind);
            cells_td = cells_td(ismember(cells_td(:,1),channels_specregion,'rows'),:);
        else
            if channels_label == brain_region
                cells_td = cells_td;
            else
                cells_td = [];
            end
        end
    end
       
    if ~isempty(cells_td)
	for icell = 1:size(cells_td,1)
        
		npage=0;
		disp(['Cell #: ' , num2str(icell)])
		if ~(dosigmoid || dorasters || doprofiles || dotuning || docellstats), break; end
		cellname = [session, num2str(cells_td(icell,1)), num2str(cells_td(icell,2))];
		
        % goodTrials_JC and goodTrials_SO numbers should be larger than 200 for each
        filename = [dirroot, cellname, '_data'];
		eval(['load ',filename])
        
        size_JC = size(goodTrials_JC,1);
        size_SO = size(goodTrials_SO,1);
        
        if size_JC < mintrialnum | size_SO < mintrialnum % 200 or 160 or 125 or 100
            disp(['Cell #: ' , num2str(icell),' fails, too few trials (at 100 trials)!'])
        else
        
            disp(['   cell ',session,' ',cellname(9:end)])
            ncells = ncells+1;
		
            if savefigs && ~onefile
                filesave = [dirroot, cellname, '_figs'];
                eval(['!del ',filesave,'.ps'])
                eval(['!del ',filesave,'.pdf'])
                close
            end
		
            %
            %sigmoid fit
            clear psyphycell
            try
                if (exnovo || savefigs), dummy; end
                filename = [dirroot, cellname, '_psyphycell'];
                eval(['load ',filename])
                [psyphycell.pairlist, psyphycell.sigmoidfit, mdl] = sigmoidfit_TT(cellname,fitt, savefigs*dosigmoid);
            catch %#ok<*CTCH>
                [psyphycell.pairlist, psyphycell.sigmoidfit, mdl] = sigmoidfit_TT(cellname,fitt, savefigs*dosigmoid);
        
                if waitshow && icell==1 && dosigmoid
                    pause(5);
                end
                
                filename = [dirroot, cellname, '_psyphycell'];
                eval(['save ',filename,' psyphycell'])
            end
                    
		
            %
            %rasters
            if dorasters && savefigs
                nplots = rasters_TT(cellname);
                for i = 1:nplots
                    figure(i); orient tall
                    eval(['print -dpsc2 -append ',filesave]) ; 
                    
                    if waitshow
                        pause(5);
                    end
				
                    if onefile
                        npage=npage+1;
                    end
				
                    close
                end
            end %if dorasters
		
            %
            %profiles
            if doprofiles
                %keyboard
                try %#ok<*UNRCH>
                    if (exnovo), dummy; end
                    filename = [dirroot, cellname, '_profiles'];
                    eval(['load ',filename])
                catch
				
                    profiles_TT(cellname);
                end
			
                if savefigs
                    nplots = profileplot_TT(cellname);
                    
                    if waitshow
                        % waitforbuttonpress
                        pause(5);
                    end                    
                    for i = 1:nplots
%                       hf=figure(i);  hf.RendererMode='manual';
%                       set(gcf,'Units','normalized', 'position',[0.5 0 0.5 1], 'PaperPositionMode','auto')
%                       eval(['print -fillpage -opengl -dpsc2 -append -r0 ',filesave]);
                        figure(i); orient tall
                        eval(['print -dpsc2 -append ',filesave]) ; 
                    
                        if waitshow
                            % waitforbuttonpress
                            pause(5);
                        end
					
                        if  onefile && i==1
                            % eval(['print -fillpage -opengl -dpsc2 -append -r0 ',filesave]);
                            npage=npage+1;
                        end
					
                        close
                    end
                end
            end %if doprofiles
		
            %
            %tuning
            if dotuning
                try
                    if (exnovo), dummy; end
                    filename = [dirroot, cellname, '_tuning'];
                    eval(['load ',filename])
                catch
                    tuning_TT(cellname);
                end
                if savefigs
                    nplots = tuningplot_TT(cellname);
                    if waitshow
						%waitforbuttonpress
                        pause(5);
                    end
                    for i = 1:nplots
%                       hf=figure(i);  hf.RendererMode='manual';
%                       set(gcf,'Units','normalized', 'position',[0.5 0 0.5 1],'Visible','off', 'PaperPositionMode','auto')
%                       eval(['print -fillpage -opengl -dpsc2 -append -r0 ',filesave]);
                        figure(i); orient tall
                        eval(['print -dpsc2 -append ',filesave]) ;
                        if waitshow
                            %waitforbuttonpress
                            pause(5);
                        end

                        if  onefile && i==1
                            % eval(['print -fillpage -opengl -dpsc2 -append -r0 ',filesave]);
                            npage=npage+1;
                        end
                        close	
                    end
                end
            end %if dotuning
		
            %cellstats
            clear cellstats anovastats tuningfit
            if docellstats
                filename = [dirroot, cellname, '_cellstats'];
                try
                    if exnovo, dummy; end
                    eval(['load ',filename])			%load cellstats file if it already exists
                    psyphycell = cellstats.psyphycell;
                    anovastats = cellstats.anovastats;
                    anovastats_both = cellstats.anovastats_both;
                    tuningfit = cellstats.tuningfit;
                    cellclass = cellstats.cellclass;
                    
                catch
                    cellstats.psyphycell = psyphycell;
                    anovastats = anovastats_TT(cellname, atleast_nntrials);
                    cellstats.anovastats = anovastats;
                    anovastats_both = anovastats_both_TT(cellname, atleast_nntrials);
                    cellstats.anovastats_both = anovastats_both;
                    cellstats.tuningfit = tuningfit_TT(cellname, anovastats, anovastats_both, atleast_nntrials, 'complete', differentRho);
                end
                filename = [dirroot,cellname,'_cellstats'];
                eval(['save ',filename,' cellstats'])
                    
                if savefigs
                    if ~checksession(cellname,.8)
                        % 						if worththestats(anovastats,.001)
                        % 						if worththestats(anovastats,.01)
                        nplots = cellstatsplot_TT(cellname, cellstats);
                        if waitshow
                            %waitforbuttonpress
                            pause(5);
                        end
                    
                        for i = 1:nplots
%                           hf=figure(i); orient tall;  hf.RendererMode='manual';
%                           set(gcf,'Units','normalized', 'position',[0.5 0 0.5 1],'Visible','off');
%                           eval(['print -opengl -dpsc2 -append -r0 ' ,filesave]);
                            figure(i); orient tall
                            eval(['print -dpsc2 -append ',filesave]) ;
                        
                            if waitshow
                                %waitforbuttonpress
                                pause(5);
                            end
                    
                            if  onefile
                                % eval(['print -besfit -opengl -dpsc2 -append -r0 ',filesave]);
                                npage=npage+1;
                            end
                            close	
                        end
                    end
                end %if savefigs
               
            end %if docellstats
	
            if mod(npage,2)
                hf=figure;
                orient tall;  hf.RendererMode='manual';
                set(gcf,'Units','normalized', 'position',[0.5 0 0.5 1],'Visible','off');
                eval(['print -opengl -dpsc2 -append -r0 ',filesave]);
            end

            close all
		
            clearvars -except npages filesaveall icell cells_td session dosigmoid dorasters doprofiles dotuning ...
			docellstats sessions dooffertest atleast_nntrials exnovo savefigs onefile fitt ...
			waitshow filesave ncells filename monkey day sessnum dirroot fileroot_ML fileroot_CED...
			filenums_CED flags cellname tic limp nfiles f filesaveall_save brain_region mintrialnum
        end %if size_JC < 200 | size_SO < 200	
        end %for icell
    end % if ~isempty(cells_td)	
end %for isession

disp(['analyzed ',num2str(ncells),' cells'])
close all

toc



