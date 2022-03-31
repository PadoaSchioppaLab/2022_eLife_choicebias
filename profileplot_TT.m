function [nplots] = profileplot_TT(cellname)
% plot activity profiles
%
% author: camillo, january 2004.
% revisions: november 2008: adapted for RC analysis
%			 december 2009: adapted for DT analysis
%			 march 2010:	adapted for new timing
%			 december 2016:	mofified for sequential offer by SB
%			 october  2018:	mofified for TT by WS

%keyboard

disp(['   ... plotting profile	of cell ',cellname])
session = cellname(1:8); readsession_TT %#ok<NASGU>
%load profiles & psyphydata
filename = [dirroot,cellname,'_profiles']; eval(['load ',filename])
trialT	= profile.specs.trial; pairs = profile.specs.pairs;
npairs = size(pairs, 1); hf = [];

%Retrieve info from profile.specs
% keyboard
allignments = profile.specs.allignments;
naligns = length(allignments);
col={'k'; 'r'; 'g'; 'b';'m';'y';'k';'c'; 'g'; 'b'; 'c';'--m';'--y';'--k';};
pairname = pairs{1}; scale=5; ms=2; zfig=0;
ncriter=numel(profile.specs.criteria); k=-naligns;

xlag=200;

selectcrit=[];
zfig=zfig+1;
hf = figure(zfig);

%
% JC average across trials
% left or right choice
crit_names = {'left','right'};
ncrits = length(crit_names);

for icrit = 1:ncrits
	crit = crit_names{icrit};
	for ialign = 1 : naligns
		eval(['prf = profile.JC.',pairname,'.side.choice.',crit,'.',allignments{ialign},';'])
		scale = max(scale,max(prf(:,2))); 
    end
end
scale=round(scale*1.5);

k=k+naligns;
for icrit = 1:ncrits
	crit = crit_names{icrit};
	for ialign = 1 : naligns
		eval(['prf = profile.JC.',pairname,'.side.choice.',crit,'.',allignments{ialign},';'])
		subplot(ncriter,naligns,k+ialign)
		hold on
		hp = plot(prf(xlag:end-xlag,1),prf(xlag:end-xlag,2),col{icrit},'Linewidth',1.5);
				
        if icrit==1
			title({ ['choice left / choice right']; [allignments{ialign}] })
			
		else
			set(gca,'ylim',[0 scale]); set(gca,'fontname','arial','fontsize',8);
			if ialign==1 && k==0
			%Cellname
				text(prf(end-100,1),scale*.95,['cell: ',cellname], 'fontsize', 8)
			elseif  ialign==naligns
				legend(crit_names{1},crit_names{2},'Offers',	'location','northeast');
			end
		end
				
		% display flagtimes
        eval(['flagtimes = profile.JC.',pairname,'.side.choice.',crit, '.flagtimes_',allignments{ialign},';'])
		for n=1:size(flagtimes,2)
			if ~isempty(flagtimes{1,n})
				hold on
				plot([flagtimes{1,n} flagtimes{1,n}],[0 scale],'k:','HandleVisibility','off')
				fill([(flagtimes{1,n}+flagtimes{2,n}) flagtimes{1,n}-flagtimes{2,n}  flagtimes{1,n}-flagtimes{2,n} (flagtimes{1,n}+flagtimes{2,n})]...
					,[scale scale 0 0],'k', 'FaceAlpha', 0.05,'linestyle','none','HandleVisibility','off');
			end
        end
        
        if strcmp(allignments{ialign},'offeron') % Add offer 1 & 2 shades
			fill([flagtimes{1,31} flagtimes{1,30}  flagtimes{1,30} flagtimes{1,31}]...
				,[scale scale 0 0],'c', 'FaceAlpha', 0.05,'linestyle','none','DisplayName','Offers');
			fill([flagtimes{1,33} flagtimes{1,32}  flagtimes{1,32} flagtimes{1,33}]...
				,[scale scale 0 0],'c', 'FaceAlpha', 0.05,'linestyle','none','HandleVisibility','off');
		end
		set(gca,'xlim',[prf(xlag,1) prf(end-xlag,1)]);
				
	end
end % end criter
        
%
% JC trial types
eval(['trialtypes = fieldnames(profile.specs.trialtypes);'])
%keyboard
%temptt=trialtypes; trialtypes(1:4)=temptt(6:9); trialtypes(5:8)=temptt(2:5); trialtypes(9)=temptt(1);
N_trialtypes = size(trialtypes,1);
for n=1:N_trialtypes
	A(n)=str2num(trialtypes{n}(2:3));
	B(n)=str2num(trialtypes{n}(5:6));
end
ordo=A./B;
[~, indtype]=sort(ordo);

zfig=zfig+1;
crit_names = {'B','A'};
ncrits = length(crit_names);
hf = figure(zfig);
z=0;
for itype = 1:N_trialtypes(1)
	trialtype = trialtypes{indtype(itype)};
	for ialign = 1:naligns
		z=z+1;
		for icrit = 1:ncrits
			crit = crit_names{icrit};
			eval(['name = profile.JC.',pairname,'.trialtype;']);
            eval(['trialtypename = ''',trialtype,'_',crit',''';']);
			if isfield(name,trialtypename)
				eval(['prf = profile.JC.',pairname,'.trialtype.',trialtype,'_',crit,'.',allignments{ialign},';'])
										
				subplot(N_trialtypes,naligns,z)
				hold on
				hp = plot(prf(xlag:end-xlag,1),prf(xlag:end-xlag,2),col{icrit},'Linewidth',1.5);
					
				% flagtime
				eval(['flagtimes = profile.JC.',pairname,'.trialtype.',trialtype,'_',crit,'.flagtimes_',allignments{ialign},';'])
			end
				
			% CONDITIONAL COSMETICS
			if icrit==1 	else
                set(gca,'ylim',[0 scale]); set(gca,'fontname','arial','fontsize',9);
                ylabel(trialtype,'fontsize',9, 'rotation',90);
                if itype==1
                    title({ ['juice choice']; [allignments{ialign}] })
				
                    if ialign==1
                    	%Cellname
                    	text(prf(end-10,1),scale*1.05,['cell: ',cellname], 'fontsize', 8)
                	elseif ialign==naligns
                    	legend(crit_names{1},crit_names{2},	'Location','best');
                    end
                end
            end
				
			% display flagtimes
			for n=1:size(flagtimes,2)
				if ~isempty(flagtimes{1,n})
					hold on
					plot([flagtimes{1,n} flagtimes{1,n}],[0 scale],'k:','HandleVisibility','off')
					fill([(flagtimes{1,n}+flagtimes{2,n}) flagtimes{1,n}-flagtimes{2,n}  flagtimes{1,n}-flagtimes{2,n} (flagtimes{1,n}+flagtimes{2,n})]...
						,[scale scale 0 0],'k', 'FaceAlpha', 0.05,'linestyle','none','HandleVisibility','off');
				end
			end
				
			if strcmp(allignments{ialign},'offeron') % Add offer 1 & 2 shades
				fill([flagtimes{1,31} flagtimes{1,30}  flagtimes{1,30} flagtimes{1,31}]...
					,[scale scale 0 0],'c', 'FaceAlpha', 0.05,'linestyle','none','DisplayName','Offers');
				fill([flagtimes{1,33} flagtimes{1,32}  flagtimes{1,32} flagtimes{1,33}]...
					,[scale scale 0 0],'c', 'FaceAlpha', 0.05,'linestyle','none','HandleVisibility','off');
			end
			set(gca,'xlim',[prf(1,1)+200 prf(end,1)-200]);
				
		end
	end
end


%
% SO average across trials
zfig=zfig+1;
hf = figure(zfig);
k=-naligns;

selectcrit=[ 5 6 4 1 2 3];
for criter=selectcrit
    crit_names = profile.specs.criteria{criter};
	ncrits = length(crit_names);
% 	SELECT BEST SCALE		
	for icrit = 1:ncrits
		crit = crit_names{icrit};
		for ialign = 1 : naligns
            eval(['prf = profile.SO.',pairname,'.',crit,'.',allignments{ialign},';'])
            scale = max(scale,max(prf(:,2))); 
		end
    end
end
scale=round(scale*1.5);

for criter=selectcrit
	k=k+naligns;
	crit_names = profile.specs.criteria{criter};
	ncrits = length(crit_names);

	for icrit = 1:ncrits
		crit = crit_names{icrit};
		for ialign = 1 : naligns
			%if criter==1 && ialign==1 % skiping meaningless results
			%	continue
			%end
				
			eval(['prf = profile.SO.',pairname,'.',crit,'.',allignments{ialign},';'])
			subplot(ncriter,naligns,k+ialign)
			hold on
			hp = plot(prf(xlag:end-xlag,1),prf(xlag:end-xlag,2),col{icrit},'Linewidth',1.5);
					
			if icrit==1
				if criter ~= 4
					title({ [crit_names{1} '/' crit_names{2}]; [allignments{ialign}] })
				else
					title({ [crit_names{1} '/' crit_names{2}]; [allignments{ialign}] })
				end
			else
				set(gca,'ylim',[0 scale]); set(gca,'fontname','arial','fontsize',8);
				if ialign==1 && k==0
					%Cellname
					text(prf(end-100,1),scale*.95,['cell: ',cellname], 'fontsize', 8)
				elseif  ialign==naligns
					legend(crit_names{1},crit_names{2},'Offers',	'location','northeast');
				end
			end
			
			% display flagtimes
			eval(['flagtimes = profile.SO.',pairname,'.',crit, '.flagtimes_',allignments{ialign},';'])
			for n=1:size(flagtimes,2)
				if ~isempty(flagtimes{1,n})
					hold on
					plot([flagtimes{1,n} flagtimes{1,n}],[0 scale],'k:','HandleVisibility','off')
					fill([(flagtimes{1,n}+flagtimes{2,n}) flagtimes{1,n}-flagtimes{2,n}  flagtimes{1,n}-flagtimes{2,n} (flagtimes{1,n}+flagtimes{2,n})]...
							,[scale scale 0 0],'k', 'FaceAlpha', 0.05,'linestyle','none','HandleVisibility','off');
				end
			end
				
			if strcmp(allignments{ialign},'offeron') % Add offer 1 & 2 shades
				fill([flagtimes{1,31} flagtimes{1,30}  flagtimes{1,30} flagtimes{1,31}]...
					,[scale scale 0 0],'c', 'FaceAlpha', 0.05,'linestyle','none','DisplayName','Offers');
				fill([flagtimes{1,33} flagtimes{1,32}  flagtimes{1,32} flagtimes{1,33}]...
					,[scale scale 0 0],'c', 'FaceAlpha', 0.05,'linestyle','none','HandleVisibility','off');
			end
			set(gca,'xlim',[prf(xlag,1) prf(end-xlag,1)]);
			
		end
	end	
end % end criter
% 	keyboard
	%autoArrangeFigures(3,1);
	
%
% SO trialtype
%profile.trialtypes.criter
eval(['trialtypes = fieldnames(profile.specs.trialtypes);'])
%keyboard
%temptt=trialtypes; trialtypes(1:4)=temptt(6:9); trialtypes(5:8)=temptt(2:5); trialtypes(9)=temptt(1);
N_trialtypes = size(trialtypes,1);
for n=1:N_trialtypes
	A(n)=str2num(trialtypes{n}(2:3));
	B(n)=str2num(trialtypes{n}(5:6));
end
ordo=A./B;
[~, indtype]=sort(ordo);

% scale=round(scale*1.2);

for criter=[1 4]
	zfig=zfig+1;
	crit_names = profile.specs.criteria{criter};
	ncrits = length(crit_names);
	hf = figure(zfig);
	z=0;
	for itype = 1:N_trialtypes(1)
		trialtype = trialtypes{indtype(itype)};
		for ialign = 1:naligns
			z=z+1;
			for icrit = 1:ncrits
				crit = crit_names{icrit};
				eval(['name = profile.SO.',pairname,'.',crit,';'])
				if isfield(name,trialtype)
					eval(['prf = profile.SO.',pairname,'.',crit,'.',trialtype,'.',allignments{ialign},';'])
					
% 					if itype==1
% 						scale = max(scale,max(prf(:,2))); scale=round(scale*1.3); % Get scale	
% 					end
					
					subplot(N_trialtypes,naligns,z)
					hold on
					hp = plot(prf(xlag:end-xlag,1),prf(xlag:end-xlag,2),col{icrit},'Linewidth',1.5);
					
					% flagtime
					eval(['flagtimes = profile.SO.',pairname,'.',crit,'.',trialtype,'.flagtimes_',allignments{ialign},';'])
				end
				
				% CONDITIONAL COSMETICS
				if icrit==1 	else
					set(gca,'ylim',[0 scale]); set(gca,'fontname','arial','fontsize',9);
					ylabel(trialtype,'fontsize',9, 'rotation',90);
					if itype==1
					if criter ~= 4
						title({ [crit_names{1} '/' crit_names{2}]; [allignments{ialign}] })
					else
						title({ [crit_names{1} '/' crit_names{2}]; [allignments{ialign}] })
					end
						if ialign==1
							%Cellname
							text(prf(end-10,1),scale*1.05,['cell: ',cellname], 'fontsize', 8)
						elseif ialign==naligns
							legend(crit_names{1},crit_names{2},	'Location','best');
						end
					end
				end
				
				% display flagtimes
				for n=1:size(flagtimes,2)
					if ~isempty(flagtimes{1,n})
						hold on
						plot([flagtimes{1,n} flagtimes{1,n}],[0 scale],'k:','HandleVisibility','off')
						fill([(flagtimes{1,n}+flagtimes{2,n}) flagtimes{1,n}-flagtimes{2,n}  flagtimes{1,n}-flagtimes{2,n} (flagtimes{1,n}+flagtimes{2,n})]...
							,[scale scale 0 0],'k', 'FaceAlpha', 0.05,'linestyle','none','HandleVisibility','off');
					end
				end
				
				if strcmp(allignments{ialign},'offeron') % Add offer 1 & 2 shades
					fill([flagtimes{1,31} flagtimes{1,30}  flagtimes{1,30} flagtimes{1,31}]...
						,[scale scale 0 0],'c', 'FaceAlpha', 0.05,'linestyle','none','DisplayName','Offers');
					fill([flagtimes{1,33} flagtimes{1,32}  flagtimes{1,32} flagtimes{1,33}]...
						,[scale scale 0 0],'c', 'FaceAlpha', 0.05,'linestyle','none','HandleVisibility','off');
				end
				set(gca,'xlim',[prf(1,1)+200 prf(end,1)-200]);
				
			end
		end
	end
end

%autoArrangeFigures(3,1);






nplots=zfig;












