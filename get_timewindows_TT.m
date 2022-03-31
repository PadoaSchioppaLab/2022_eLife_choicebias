function [timewindows, colrs, timewindows_off, timewindows_cho, timewindows_tgt, timewindows_cue,  timewindows_bas] = get_timewindows_TT(JCorSO)
%
% defining time windows
%
if isequal(JCorSO, 'JC')
    % % win number		win name		ref flags	ref times
    timewindows{1}	=	{'preoffer',	[30, 30],	[-500,   0]};
    timewindows{2}	=	{'postoffer',	[30, 30],	[ 100, 600]};
    timewindows{3}	=	{'latedelay',	[30, 30],	[ 600,1100]};
    timewindows{4}	=	{'prego',		[35, 35],	[-400, 100]};
    timewindows{5}	=	{'reactime',	[35, 39],	[   0,   0]};
    timewindows{6}	=	{'prejuice',	[50, 50],	[-400, 100]};
    timewindows{7}	=	{'postjuice',	[50, 50],	[ 100, 600]};
    
    % super detailed twds
    start=-500; step=25;
    for n=1:181
        timewindows_off{n}	=	{['twds_off' num2str(n)],	[30, 30],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end

    % super detailed twds & OUTCOME
    start=-1500; step=25;
    for n=1:181
        timewindows_cho{n}	=	{['twds_cho' num2str(n)],	[50, 50],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end

    % super detailed twds % TGT
    start=-1500; step=25;
    for n=1:121
        timewindows_tgt{n}	=	{['twds_tgt' num2str(n)],	[35, 35],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end
    
    % super detailed twds % CUE
    start=-100; step=25;
    for n=1:41
        timewindows_cue{n}	=	{['twds_cue' num2str(n)],	[25, 25],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end
    
    % super detailed twds % baseline(BAS)
    start=-500; step=25;
    for n=1:1
        timewindows_bas{n}	=	{['twds_bas' num2str(n)],	[25, 25],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end
    
    colrs = [
        .75 .75 .75;%grey		posttrialcue
        0   0  .7;	%dark blue	postoffer
        1   0   1;	%magenta	latedelay
        1   0   0;	%red		prego
        0	.7	0;	%green		reactime
        0 .75 .75;	%celeste	prejuice
        0   0   0	%black		postjuice
        ];	

elseif isequal(JCorSO, 'SO')
    
    %500ms + late postjuice + physiological delays
    % win number		win name		ref flags	ref times
    timewindows{1}	=	{'preoffers',	[30, 30],	[-500,   0]};
    timewindows{2}	=	{'postoffer1',	[30, 30],	[   100, 600]};
    timewindows{3}	=	{'interoffers',	[32, 32],	[ -400, 100]};
    timewindows{4}	=	{'postoffer2',	[32, 32],	[  100, 600]};
    timewindows{5}	=	{'offersoff',	[32, 32],	[	600, 1100]};
    timewindows{6}	=	{'fixoff',		[27, 27],	[	-100, 400]};
    timewindows{7}	=	{'prejuice',	[50, 50],	[	-400, 100 ]};
    timewindows{8}	=	{'postjuice',	[50, 50],	[  100, 600]};
    timewindows{9}	=	{'postjuice2',	[50, 50],	[	600, 1100]}; %% 

    % super detailed twds & OFFER 1 OR OFFER 2
    % OFFER 1
    start=-500; step=25;
    for n=1:181
        timewindows_off{n}	=	{['twds_off' num2str(n)],	[30, 30],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end

    % super detailed twds & OUTCOME
    start=-1500; step=25;
    for n=1:181
        timewindows_cho{n}	=	{['twds_cho' num2str(n)],	[50, 50],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end

    % super detailed twds % TGT
    start=-1500; step=25;
    for n=1:121
        timewindows_tgt{n}	=	{['twds_tgt' num2str(n)],	[27, 27],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end
 
    % super detailed twds % CUE
    start=-100; step=25;
    for n=1:41
        timewindows_cue{n}	=	{['twds_cue' num2str(n)],	[25, 25],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end

    % super detailed twds % baseline(BAS)
    start=-500; step=25;
    for n=1:1
        timewindows_bas{n}	=	{['twds_bas' num2str(n)],	[25, 25],	[start+((n-1)*step), start+((n-1)*step)+100]};
    end
    
    colrs = [
        .6 .6 .6;	%grey		preoffers/posttrialcue
        .6   0  .6;	%light magenta	offer1
        0   0  .5;	%light blue	interoffers
        1   0   1;	%magenta	offer2
        0   0  .8;	%blue		offersoff
        0	.8	0;	%green		reactime
        0	.8 .8;	%celeste	prejuice
        .5	0	0;	%light red	postjuice
        .8	0	0;	%red		postjuice2
        ];
    
end

