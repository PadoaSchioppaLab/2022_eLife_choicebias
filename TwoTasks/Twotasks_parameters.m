
function [sessionParams] = Twotasks_parameters()
%
% parameters for SequentialOffers task
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Mode parameters %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sessionParams.fakemonkey	= 0;
sessionParams.p_nofixbreak	= 1;	%.99;
sessionParams.manualchoice	= 0;	%press Home/PgUp key for L/R choice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% State system parameters %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sessionParams.maxNoTrials = 2000;		%max number of trials
sessionParams.maxNoBlocks = 1000;			%max number of blocks

%define time intervals (in ms):
% SO
sessionParams.trial.wait_for_fix	= 4000; %5000
sessionParams.trial.initial_prefix	= 500;
sessionParams.trial.initial_fix_SO  = 500; %950
sessionParams.trial.initial_fix_JC  = 500; %1000
sessionParams.trial.initial_idle	= 200;
sessionParams.trial.offer1on_time	= 500; %650
sessionParams.trial.interoffer_time	= 500;
sessionParams.trial.offer2on_time	= 500; %650
sessionParams.trial.offeroff_time	= 500;
sessionParams.trial.tgton_range     = [500 1000];
sessionParams.trial.toofast         = 10;
sessionParams.trial.max_reactime	= 1000;
sessionParams.trial.max_sacc_time	= 1000;
sessionParams.trial.sacctgt_idle	= 200;
sessionParams.trial.tgthold_time	= 600;
sessionParams.trial.preend_time     = 1000;

%JC
sessionParams.trial.delay_range		= [1000 2000];
% sessionParams.trial.tgton_range_JC  = [1500 2000];
sessionParams.stimulus.offerDist	= 7;	%distance from center for offer stimuli

%fixation window (in degrees):
sessionParams.prefix_radius	= 4;
sessionParams.fix_radius	= 6; %4
%sessionParams.tgt_radius	= 7; %5
sessionParams.tgt_radiusL	= 4; %5
sessionParams.tgt_radiusR	= 5; %5


sessionParams.tries = 10;

sessionParams.rand_outcome.psudorand = 1;
sessionParams.rand_outcome.phat0 = .10;		%prob threshold for poor luck
sessionParams.rand_outcome.phat1 = .05;		%prob threshold for good luck

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Parameters for goods %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% goods
% goods(0).juice  		= 'lemon k.a.';
% goods(0).juice  		= 'grape';
% goods(0).juice  		= '3/4 peach';
% goods(0).juice  		= '1/3 cranberry';
% goods(0).juice  		= '2/3 fruit punch';
% goods(0).juice  		= '1/2 apple';
% goods(0).juice  		= 'cranberry';
% goods(0).juice  		= 'lemon k.a.';
% goods(0).juice  		= 'peppermint';
% goods(0).juice  		= 'water';
% goods(0).juice  		= 'Tamarind k.a.';
% goods(0).juice  		= '.65 g/l salt';
% goods(0).juice  		= 'watermelon';
% goods(0).juice  		= 'h punch';
% goods(0).juice  		= 'cherry';

%

%
goods(1).juice  		= 'lemon k.a.';
goods(1).juiceline  	= 3;

goods(2).juice  		= 'grape';
goods(2).juiceline  	= 1;


%
ngoods = size(goods,2);
for i = 1:ngoods
    if ~isempty(goods(i).juice)
        goods(i).probability= 1;	%#ok<AGROW>
        goods(i).quantum	= 71;	%#ok<AGROW>
        goods(i).cost		= 0;	%#ok<AGROW>
        %
        if		isequal(goods(i).juice, 'lemon k.a.'),		goods(i).color = [1 1 0];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, '3/4 peach'),		goods(i).color = [1 .4 .4];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, 'grape'),			goods(i).color = [0 1 0];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, '1/2 apple'),		goods(i).color = [1 .66 0];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, 'berry blue'),		goods(i).color = [.5 .5 1];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, '3/4 cherry'),			goods(i).color = [.85 0 0];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, '2/3 fruit punch'),	goods(i).color = [1 0 .65];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, 'h punch'),			goods(i).color = [.5 .5 1];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, '1/3 cranberry'),   goods(i).color = [1 .7 .85];	%#ok<AGROW>
        elseif	isequal(goods(i).juice, 'peppermint'),		goods(i).color = [0 .8 .8];		%#ok<AGROW>
            % 		elseif	isequal(goods(i).juice, 'water'),			goods(i).color = [1 1 1];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, 'water'),			goods(i).color = [0 1 1];		%#ok<AGROW>
            % 		elseif	isequal(goods(i).juice, 'watermellon'),		goods(i).color = [1 .2 .3];		%#ok<AGROW>
        elseif	isequal(goods(i).juice, 'watermelon'),		goods(i).color = [.55 .8 .55];	%#ok<AGROW>
        elseif	isequal(goods(i).juice, 'Tamarind k.a.'),	goods(i).color = [.71 .53 .53];	%#ok<AGROW>
        elseif	isequal(goods(i).juice, '.65 g/l salt'),	goods(i).color = [.75 .75 .75];	%#ok<AGROW>
        else		disp('good(',num2str(i),').color not defined!!'), dummy
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Symbols representing the temporal order %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% symbols.offer1 = 'sqr';
% symbols.offer2 = 'crc';

% symbols.offer1 = 'square';
% symbols.offer2 = 'circle';

symbols.offer1 = 'square';
symbols.offer2 = 'square';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Offer list %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%TT%%%%
%1: JC
%2: SO

% % %
% offerTypes = [  %1A~=2.5B %for training
% %	%N	#A	#B
%     2   3   0   1
%     2   4   0   1
%     2   1   5   1
%     2   2   4   1
%     2   2   8   1
%     2   3   6   1
%     2   3  10   1
%     2   4  10   1
%     2   0   5   1
%     2   0   6   1
%     
%     2   3   0   2
%     2   4   0   2
%     2   1   5   2
%     2   2   4   2
%     2   2   8   2
%     2   3   6   2
%     2   3  10   2
%     2   4  10   2
%     2   0   5   2
%     2   0   6   2
% 	];

% % % 
% offerTypes = [	%1A~=1.5B; A:0-4; B:0-6; standard edited for SO
% %	%N	#A	#B
% 	1	3	0   1
% 	2	1	2   1
% 	2	1	4   1
% 	2	2	2   1
% 	2	2	4   1
% 	2	2	6   1
% 	2	3	4   1
%   2   3   6   1
%   2   4   1   1
% 	2	4	3   1
%   2   4   6   1
% 	1	0	4   1
%     
%   1	3	0   2
% 	3	1	2   2
% 	2	1	4   2
% 	3	2	2   2
% 	3	2	4   2
% 	3	2	6   2
% 	3	3	4   2
%   3   3   6   2
%   2   4   1   2
% 	3	4	3   2
%   3   4   6   2
%   1	0	4   2
% 	];

% % % 
% offerTypes = [	%1A~=1.5B; A:0-4; B:0-6; standard edited for SO; OVA&OVB anti-correlated
% %	%N	#A	#B
% 	1	3	0   1
% 	2	1	4   1
%     2   1   6   1
% 	2	2	2   1
%     2   2   3   1
% 	2	2	4   1
%     2   2   6   1
%     2   3   1   1
%     2   3   3   1
% 	2	3	4   1
%     2   4   1   1
% 	2	4	3   1
% 	1	0	4   1
%     
%     1	3	0   2
% 	2	1	4   2
%     2   1   6   2
% 	2	2	2   2
%     2   2   3   2
% 	2	2	4   2
%     2   2   6   2
%     2   3   1   2
%     2   3   3   2
% 	2	3	4   2
%     2   4   1   2
% 	2	4	3   2
% 	1	0	4   2
% 	];

% % % 
% offerTypes = [	%1A~=1.5B; A:0-6; B:0-8; standard edited for SO
% %	%N	#A	#B
% 	1	3	0   1
% 	2	2	2   1
% 	2	2	4   1
% 	2	2	6   1
%   2   3   8   1
%   2   4   4   1
% 	2	4	8   1
%   2   6   2   1
%   2   6   6   1
%   2   6   8   1
% 	1	0	4   1
%     
%   1	3	0   2
% 	3	2	2   2
% 	3	2	4   2
% 	3	2	6   2
%   3   3   8   2
%   3   4   4   2
% 	3	4	8   2
%   3   6   2   2
%   3   6   6   2
%   3   6   8   2
% 	1	0	4   2
% 	];

% % % 
offerTypes = [  %1A~=1.5B; A:0-6; B:0-8; standard edited for SO; OVA&OVB anti-correlated
%	%N	#A	#B
    1	3   0   1
    2   1   8   1
	2	2	4   1
    2   2   6   1
    2   2   8   1
    2   3   3   1
    2   3   4   1
    2   3   6   1
    2   4   2   1
    2   4   3   1
    2   4   4   1
    2   4   6   1
    2   6   1   1
    2   6   2   1
    % 2   6   8   1
	1   0   4   1
    
    1	3	0   2
    2   1   8   2
	2	2	4   2
    2   2   6   2
    2   2   8   2
    2   3   3   2
    2   3   4   2
    2   3   6   2
    2   4   2   2
    2   4   3   2
    2   4   4   2
    2   4   6   2
    2   6   1   2
    2   6   2   2
    % 2   6   8   2
	1	0	4   2
	];

% % %
% offerTypes = [	%1A~=2B; A:0-4; B:0-10; standard edited for SO
% %	%N	#A	#B  TT
%     1   3   0   1   
%     2   1   8   1
%     2   2   2   1
%     2   2   6   1
%     2   3   4   1
%     2   3  10   1
%     2   4   2   1
%     2   4   6   1
%     2   4   8   1
%     2   4  10   1
%     1   0   4   1
%     
%     1   3   0   2   
%     3   1   8   2
%     3   2   2   2
%     3   2   6   2
%     3   3   4   2
%     3   3  10   2
%     3   4   2   2
%     3   4   6   2
%     3   4   8   2
%     3   4  10   2
%     1   0   4   2
% 	];

% % %
% offerTypes = [	%1A~=2B; A:0-4; B:0-10; standard edited for SO; OVA&OVB anti-correlated
% %	%N	#A	#B  TT
%     1   3   0   1
%     2   1   4   1
%     2   1   8   1
%     2   2   4   1
%     2   2   6   1
%     2   2   8   1
%     2   2   10  1
%     2   3   2   1
%     2   3   4   1
%     2   3   6   1
%     2   4   2   1
%     2   4   4   1
%     1   0   4   1
%     
%     1   3   0   2
%     2   1   4   2
%     2   1   8   2
%     2   2   4   2
%     2   2   6   2
%     2   2   8   2
%     2   2   10  2
%     2   3   2   2
%     2   3   4   2
%     2   3   6   2
%     2   4   2   2
%     2   4   4   2
%     1   0   4   2
% 	];

% % % 
% offerTypes = [	%1A~=2B; A:0-3; B:0-8; standard edited for SO
% %	%N	#A	#B
% 	1	3	0   1
%     2	1	3   1
%     2	1	8   1
%     2	2	3   1
% 	2	2	4   1
% 	2	2	6   1
%     2   3   4   1
%     2	3	6   1
% 	2	3	8   1
% 	1	0	4   1
% 
% 	1	3	0   2
%     3	1	3   2
%     3	1	8   2
%     3	2	3   2
% 	3	2	4   2
% 	3	2	6   2
%     3   3   4   2
%     3	3	6   2
% 	3	3	8   2
% 	1	0	4   2
% 	];

% % % 
% offerTypes = [	%1A~=2B; A:0-3; B:0-8; standard edited for SO; OVA&OVB anti-correlated
% %	%N	#A	#B
% 	1	3	0   1
%     2	1	6   1
%     2	1	8   1
%     2   2   2   1
%     2	2	3   1
% 	2	2	4   1
% 	2	2	6   1
%     2   3   2   1
%     2   3   3   1
%     2   3   4   1
%     2	3	6   1
% 	1	0	4   1
% 
% 	1	3	0   2
%     3	1	6   2
%     3	1	8   2
%     3   2   2   2
%     3	2	3   2
% 	3	2	4   2
% 	3	2	6   2
%     3   3   2   2
%     3   3   3   2
%     3   3   4   2
%     3	3	6   2
% 	1	0	4   2
% 	];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aux = sum(offerTypes,1)>0;
jj = find(aux,1,'last');
offerTypes = offerTypes(:,1:jj);

%
% condense offType matrix and check that in any trial noff <= 2
ind = offerTypes(:,1)>0;
offtypes = offerTypes(ind,:);
aux = sum(offtypes(:,2:end-1)>0,2);
if (max(aux)>2), disp('problem: offerType with >2 goods!!'); dummy; end
offerList = [];
for itype = 1:size(offtypes,1)
    for iiN = 1:offtypes(itype,1)
        offerList = [offerList; offtypes(itype,2:end)];		%#ok<AGROW>
    end
end

ngoods = size(offerList,2)-1;
goods = goods(1:ngoods);

% for the records
sessionParams.goods			= goods;
sessionParams.symbols		= symbols;
sessionParams.offerTypes	= offerTypes;
sessionParams.offerList		= offerList;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Parameters for visual stimuli %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sessionParams.stimulus.innerRadius	= 2.0;	%inner radius for stimuli
% sessionParams.stimulus.outerRadius	= 4.0;	%outer radius for stimuli
sessionParams.stimulus.innerRadius	= 2.5;	%inner radius for stimuli
sessionParams.stimulus.outerRadius	= 5.0;	%outer radius for stimuli
% sessionParams.stimulus.offerDist	= 7.0;	%distance from center for offer stimuli
sessionParams.stimulus.squareSize	= 1.0;	%size in deg
sessionParams.stimulus.diamondSize	= 1.0;
sessionParams.stimulus.crossSize	= 1.0;
sessionParams.stimulus.circleSize	= 1.0;
sessionParams.stimulus.lineSize		= 1.0;
% sessionParams.stimulus.targetDist	= 7.0;	%distance from center for saccade targets
sessionParams.stimulus.targetDist	= 7.0;	%distance from center for saccade targets
sessionParams.stimulus.targetSize	= .25;
sessionParams.stimulus.targetFill	= 1;	%0 if not filled
sessionParams.stimulus.targetColor	= [1 1 1];
sessionParams.stimulus.nTgtPositions= 2;


%
symbPos_all_JC	= []; symbPos_all_SO	= [];
inRadius		= sessionParams.stimulus.innerRadius;
outRadius		= sessionParams.stimulus.outerRadius;
for i = 1:6,	symbPos_all_JC = [symbPos_all_JC; inRadius*[cos((i-1)*pi/3),sin((i-1)*pi/3)]]; end		%#ok<AGROW>
for i = 1:12,	symbPos_all_JC = [symbPos_all_JC; outRadius*[cos((i-1)*pi/6),sin((i-1)*pi/6)]]; end		%#ok<AGROW>
for i = 1:6,	symbPos_all_SO = [symbPos_all_SO; inRadius*[cos((i-1)*pi/3 +pi/2),sin((i-1)*pi/3 +pi/2)]]; end		%#ok<AGROW>
for i = 1:12,	symbPos_all_SO = [symbPos_all_SO; outRadius*[cos((i-1)*pi/6 +pi/2),sin((i-1)*pi/6 +pi/2)]]; end		%#ok<AGROW>
% %
stimIndices{1} = 1;
stimIndices{2} = [2 6]';
stimIndices{3} = [1 2 6]';
stimIndices{4} = [2 3 5 6]';
stimIndices{5} = [1 2 6 10 16]';
stimIndices{6} = [2 3 5 6 10 16]';
stimIndices{7} = [1 2 3 5 6 10 16]';
stimIndices{8} = [2 3 5 6 9 11 15 17]';
stimIndices{9} = [2 3 5 6 9 10 11 15 17]';
stimIndices{10}= [2 3 5 6 9 10 11 15 16 17]';
stimIndices{11}= [1 2 3 5 6  9 10 11 15 16 17]';
stimIndices{12}= [2 3 5 6 8  9 11 12 14 15 17 18]';
stimIndices{13}= [1 2 3 5 6  8  9 11 12 14 15 17 18]';
stimIndices{14}= [2 3 5 6 8  9 10 11 12 14 15 16 17 18]';
stimIndices{15}= [1 2 3 5 6  8  9 10 11 12 14 15 16 17 18]';
stimIndices{16}= [1 2 3 5 6  7  8  9 10 11 12 14 15 16 17 18]';
%
sessionParams.symbPos_all_JC = symbPos_all_JC;
sessionParams.symbPos_all_SO = symbPos_all_SO;
sessionParams.stimIndices = stimIndices;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Other parameters %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %
% % sounds
% samplefreq = 8192;
% sessionParams.sounds.samplefreq = samplefreq;
% %
% % poorlucksound
% sound_ms = 250;
% xy = 100;
% seg = pi*(1:xy)'/xy;
% ii = 8;
% wf = [];
% for jj = 1 : sound_ms/1000 * samplefreq/xy
% 	wf = [wf; sin(ii*seg)];		%#ok<AGROW>
% end
% sessionParams.sounds.poorluck = wf;
% %sound(wf,samplefreq)

