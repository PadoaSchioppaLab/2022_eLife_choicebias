function [C timingfile TrialRecord userdef_trialholder_function] = Twotasks_tlf(MLConfig, TrialRecord)
%
%
%
userdef_trialholder_function = [];

%persistent sessionParams offerListHalfBlock
persistent TwoID offerListHalfBlock targetListHalfBlock		%%to fix: should be saved??

%timing file
timingfile = 'Twotasks_timingfile.m';

%load parameters used in this session
sessionParams = Twotasks_parameters;
offerList = sessionParams.offerList;
nTgtPositions = sessionParams.stimulus.nTgtPositions;
ngoods = size(offerList,2);

%
BackgroundColor=MLConfig.SubjectScreenBackground;
DAQ = MLConfig.DAQ;
ScreenInfo = MLConfig.Screen;
% keyboard
trialNumber = TrialRecord.CurrentTrialNumber;
testtrial = 0;

if ~trialNumber %initialization trial
    trialNumber = 1;
    testtrial = 1;
    % 	juiceDeliver(DaqInfo, NaN, 0, 0, -1);
end

if trialNumber==1
    TrialRecord.User.sessionParams				= sessionParams;
    TrialRecord.User.ngoods						= ngoods;
    TrialRecord.User.blocksDone					= 0;
    TrialRecord.User.halfBlocksDone				= 0;
    TrialRecord.User.trialsDoneThisHalfBlock		= 0;
    TrialRecord.User.goodTrials_JC					= [];
    TrialRecord.User.goodTrials_SO					= [];
    TrialRecord.User.summaryStats.choice12		= [0 0];
    TrialRecord.User.summaryStats.choiceLR		= [0 0];
    TrialRecord.User.summaryStats.gotjuice		= zeros(2,ngoods);
    TrialRecord.User.summaryStats.gotjuice_hist	= NaN(sessionParams.maxNoTrials,ngoods);
    TrialRecord.User.summaryStats.psudorand		= zeros(sessionParams.maxNoTrials,1);
    TrialRecord.User.choicePattern				= [];
    
    TrialRecord.User.summaryStats.LRchoice		= [0 0];
    TrialRecord.User.summaryStats.quantity		= zeros(1,ngoods);
    
    
end
%
blocksDone				= TrialRecord.User.blocksDone;
halfBlocksDone			= TrialRecord.User.halfBlocksDone;
trialsDoneThisHalfBlock	= TrialRecord.User.trialsDoneThisHalfBlock;

%
%if we're done, return
if		isequal(blocksDone, sessionParams.maxNoBlocks)
    disp('... completed max number of blocks ...')
    return
elseif	trialNumber > sessionParams.maxNoTrials
    disp('... completed max number of trials ...')
    return
end

%if beginning of half-block, generate offerListHalfBlock
if trialNumber == 1 || halfBlocksDone == 2		% block ended
    
    offerListHalfBlock = [];
    for itgtpos = 1:nTgtPositions
        offerListHalfBlock = [offerListHalfBlock; [offerList, itgtpos*ones(size(offerList,1),1)]]; %#ok<AGROW>
    end
    % targetListHalfBlock = offerListHalfBlock(:,3);
    targetListHalfBlock = offerListHalfBlock(:,4); % WS
    [offtypes, junk, offtypenum] = unique(offerListHalfBlock(:,[1 2 3]),'rows');
    TwoID = offerListHalfBlock(:,3);
%      keyboard
    
    orderA = ones(size(offerListHalfBlock,1),1);	%order of offer A
    for iofftype = 1:size(offtypes,1)
        ind = find(offtypenum==iofftype);
        ind = ind(randperm(length(ind)));
        ind = ind(1:round(length(ind)/2));
%           if TwoID(iofftype)==2 % DO IT ONLY FOR SEQ OFFER
        orderA(ind) = -orderA(ind);				%here the sign represent the temporal order: -/+ for offer 1/2
%           end
     end
%     TwoID = offerListHalfBlock(:,3);
    offerListHalfBlock = offerListHalfBlock(:,1:2).*[orderA,-orderA];
    nTrialsHalfBlock = size(offerListHalfBlock,1);
    
    halfBlocksDone = 0;
    trialsDoneThisHalfBlock = zeros(nTrialsHalfBlock,1);
end
nTrialsHalfBlock = size(offerListHalfBlock,1);

%if halfBlock just ended
if sum(trialsDoneThisHalfBlock) == nTrialsHalfBlock
    trialsDoneThisHalfBlock = zeros(nTrialsHalfBlock,1);
    offerListHalfBlock = -offerListHalfBlock;
end
%
blockNumber = blocksDone + 1;

%
%disp some info
choice12 = TrialRecord.User.summaryStats.choice12;
choiceLR = TrialRecord.User.summaryStats.choiceLR; % SO
LRchoice = TrialRecord.User.summaryStats.LRchoice; % JC
nhits = sum(choiceLR); nhits2 = sum(LRchoice);
perc12 = round(100*choice12/nhits);
percLR = round(100*choiceLR/nhits);
LRperc = round(100*LRchoice/nhits2);
%
disp(['trial : ',num2str(nhits),' of ',num2str(trialNumber)])
disp(['  blocks done: ',num2str(blocksDone), '  halfblocks done: ',num2str(halfBlocksDone), ...
    '  trials done this halfblock: ',num2str(sum(trialsDoneThisHalfBlock)),' of ',num2str(nTrialsHalfBlock)])
disp(['%choices 1/2 (SO): ', num2str(perc12)])
disp(['%choices L/R (SO): ', num2str(percLR)])
disp(['%choices L/R (JC): ', num2str(LRperc)])
%
jrec = TrialRecord.User.summaryStats.gotjuice;
str = '';
for ig = 1:ngoods
    str = [str,num2str(jrec(2,ig)),'/',num2str(jrec(1,ig)),'  ']; %#ok<AGROW>
end
disp(['    got juice: ',str])
disp('  ')

%initialize level vars to zero
% keyboard
putvalue(DAQ.TTL{1}, 0);
% putvalue(DaqInfo.TTL2, 0);
%[offerListHalfBlock]
%dummy
%
% randomly select next offer type (without replacement)
posstrials = find(trialsDoneThisHalfBlock==0);
ind = ceil(rand(1)*length(posstrials));
trialIndex = posstrials(ind);
currentOfferType = offerListHalfBlock(trialIndex,:);
currentTgtNumber = targetListHalfBlock(trialIndex);
%
otherTgtNumber = mod(currentTgtNumber + nTgtPositions/2, nTgtPositions);
if ~otherTgtNumber, otherTgtNumber = nTgtPositions; end
%
% compute target position associated to offer(1) (left)
if nTgtPositions==2
    th = 2*pi*(currentTgtNumber-1)/nTgtPositions;	%puts the 2 targets left/right -- remove for up/down
else
    th = 2*pi*currentTgtNumber/nTgtPositions - pi/nTgtPositions;
end
xytargetposition = double(sessionParams.stimulus.targetDist*[cos(th), sin(th)]);
%xytargetposition = [-sessionParams.stimulus.offerDist, 0];	% by the offer

%
% current offers
currentOffers = [];

%  TwoID(trialIndex)==1 % 1=SO 2=JC
currentOffers(1).Twotasks = TwoID(trialIndex);
currentOffers(2).Twotasks = TwoID(trialIndex);

ind1 = find(currentOfferType<0);	%good associated with offer 1
ind2 = find(currentOfferType>0);	%good associated with offer 2

if ~isempty(ind1)
    currentOffers(1).goodId			= ind1;
    currentOffers(1).good			= sessionParams.goods(ind1);
    currentOffers(1).quantity		= abs(currentOfferType(ind1));
else
    currentOffers(1).goodId			= 0;
    currentOffers(1).good.color		= sessionParams.goods(setdiff([1 2],ind2)).color;
    currentOffers(1).quantity		= 0;
end
currentOffers(1).order				= 1;
currentOffers(1).symbol				= sessionParams.symbols.offer1;
if currentOffers(1).Twotasks==1
    currentOffers(1).xyposition		= [-sessionParams.stimulus.offerDist, 0];
elseif currentOffers(1).Twotasks==2
    currentOffers(1).xyposition			= [0, 0];		%all offers are centered
end

if currentOffers(1).Twotasks==1
    % currentOffers(1).xytargetposition	= -xytargetposition;
    % currentOffers(1).targetnumber		= currentTgtNumber;
    currentOffers(1).targetnumber		= 1;
    currentOffers(1).xytargetposition	= currentOffers(1).xyposition;
elseif currentOffers(1).Twotasks==2
    currentOffers(1).targetnumber		= currentTgtNumber;
    currentOffers(1).xytargetposition	= -xytargetposition;
end

if ~isempty(ind2)
    currentOffers(2).goodId			= ind2;
    currentOffers(2).good			= sessionParams.goods(ind2);
    currentOffers(2).quantity		= currentOfferType(ind2);
else
    currentOffers(2).goodId			= 0;
    currentOffers(2).good.color		= sessionParams.goods(setdiff([1 2],ind1)).color;
    currentOffers(2).quantity		= 0;
end
currentOffers(2).order				= 2;
currentOffers(2).symbol				= sessionParams.symbols.offer2;
if currentOffers(2).Twotasks==1
    currentOffers(2).xyposition		= [+sessionParams.stimulus.offerDist, 0];
elseif currentOffers(2).Twotasks==2
    currentOffers(2).xyposition			= [0, 0];		%all offers are centered
end

% currentOffers(2).targetnumber		= otherTgtNumber;
if currentOffers(2).Twotasks==1
    % currentOffers(1).xytargetposition	= -xytargetposition;
    % currentOffers(2).targetnumber		= otherTgtNumber;
    currentOffers(2).targetnumber		= 2;
    currentOffers(2).xytargetposition	= currentOffers(2).xyposition;
elseif currentOffers(2).Twotasks==2
    currentOffers(2).targetnumber		= otherTgtNumber;
    currentOffers(2).xytargetposition	= +xytargetposition;
end

prefix_point    = 1;
fix_point_SO	= 2;
fix_point_JC	= 3;
fix_point_dim_JC= 4;
fixpoint_offer1	= 5;
fixpoint_offer2	= 6;
target1         = 7;
target2     	= 8;
offer1      	= 9;
offer2      	= 10;


%
% generate task objects for this trial
C = [];


% prefix point
C(1).Type	= 'crc';
C(1).Name	= 'prefixationPoint';
C(1).Xpos	= 0;
C(1).Ypos	= 0;
C(1).Radius = .5;
C(1).Color = [1 1 1];
C(1).FillFlag = 1;

% fix point_SO
C(2).Type	= 'crc';
C(2).Name	= 'fixationPoint_SO';
C(2).Xpos	= 0;
C(2).Ypos	= 0;
C(2).Radius = .2;
C(2).Color = [1 1 1];
C(2).FillFlag = 1;

% fix point_JC
fix1= (ones(111,111,3)); %%% CREATE CROSS
for n=1:3; fix1(:,:,n)=fix1(:,:,n).*BackgroundColor(n); end;
fix1(40:70,:,:)=1; fix1(:,40:70,:)=1;
imwrite(fix1,'fixationPoint_JC','jpeg');
C(3).Type	= 'pic';
C(3).Name	= 'fixationPoint_JC';
C(3).Xpos	= 0 ; %centerLocation(1);
C(3).Ypos	= 0 ; %centerLocation(2);
C(3).Xsize = 20;
C(3).Ysize = 20;


% % dim fix point (BackgroundColor)
C(4) = C(3);
fix1= (ones(111,111,3)); %%% CREATE CROSS
for n=1:3; fix1(:,:,n)=fix1(:,:,n).*BackgroundColor(n); end;
fix1(40:70,:,:)=BackgroundColor(n); fix1(:,40:70,:)=BackgroundColor(n);
imwrite(fix1,'fixationPoint_dim','jpeg');
C(4).Name	= 'fixationPoint_dim';
% C(4).Color = ([1 1 1] + BackgroundColor)/2;


% color donuts around fix point
C(5)=C(2);
C(5).Name	= 'fixationPoint_Offer1';
C(5).Xpos	= 0;
C(5).Ypos	= 0;
C(5).Radius = .55;
C(5).Color = [currentOffers(1).good.color];
C(5).FillFlag = 0;

% % color donuts around fix point
C(6)=C(2);
C(6).Name	= 'fixationPoint_Offer2';
C(6).Xpos	= 0;
C(6).Ypos	= 0;
C(6).Radius = .55;
C(6).Color = [currentOffers(2).good.color];
C(6).FillFlag = 0;

warning off MATLAB:intConvertNonIntVal
%
% saccade targets
C(7).Type	= 'crc';
C(7).Name	= 'targetImage1';
C(7).Xpos	= currentOffers(1).xytargetposition(1);
C(7).Ypos	= currentOffers(1).xytargetposition(2);
C(7).Radius = sessionParams.stimulus.targetSize;
if currentOffers(1).Twotasks==1
    C(7).Color	= [1 1 1];
elseif currentOffers(1).Twotasks==2
    C(7).Color	= currentOffers(1).good.color;
end
C(7).FillFlag = sessionParams.stimulus.targetFill;

%
C(8).Type	= 'crc';
C(8).Name	= 'targetImage2';
C(8).Xpos	= currentOffers(2).xytargetposition(1);
C(8).Ypos	= currentOffers(2).xytargetposition(2);
C(8).Radius = sessionParams.stimulus.targetSize;
if currentOffers(2).Twotasks==1
    C(8).Color	= [1 1 1];
elseif currentOffers(2).Twotasks==2
    C(8).Color	= currentOffers(2).good.color;
end
C(8).FillFlag = sessionParams.stimulus.targetFill;

%
% visual stimuli
if currentOffers(1).Twotasks==1
    offerImage = juicechoice_makestimuli(currentOffers(1), sessionParams, MLConfig);
elseif currentOffers(1).Twotasks==2
    offerImage = SequentialOffers_makestimuli(currentOffers(1), sessionParams, MLConfig, 'offer');
end
imwrite(offerImage,'offerImage1','jpeg');
C(9).Type	= 'pic';
C(9).Name	= 'offerImage1';
C(9).Xpos	= currentOffers(1).xyposition(1);
C(9).Ypos	= currentOffers(1).xyposition(2);

%
if currentOffers(2).Twotasks==1
    offerImage = juicechoice_makestimuli(currentOffers(2), sessionParams, MLConfig);
elseif currentOffers(2).Twotasks==2
    offerImage = SequentialOffers_makestimuli(currentOffers(2), sessionParams, MLConfig, 'offer');
end
imwrite(offerImage,'offerImage2','jpeg');
C(10).Type	= 'pic';
C(10).Name	= 'offerImage2';
C(10).Xpos	= currentOffers(2).xyposition(1);
C(10).Ypos	= currentOffers(2).xyposition(2);
warning on MATLAB:intConvertNonIntVal


for i=1:length(C)
    if strcmp(C(i).Type,'pic') && ~strcmp(C(i).Name,'fixationPoint_JC') && ~strcmp(C(i).Name,'fixationPoint_dim');
        C(i).Xsize = size(offerImage,1); % FIXSB
        C(i).Ysize = size(offerImage,2); % FIXSB
        %     elseif ~strcmp(C(i).Name,'fixationPoint_JC') || ~strcmp(C(i).Name,'fixationPoint_dim') ;
        %         C(i).Xsize = -1;
        %         C(i).Ysize = -1;
    end
end

%
% update records
currTrial.Twotasks	= [];
currTrial.trialNumber	= trialNumber;
currTrial.blockNumber	= blockNumber;
currTrial.trialIndex	= trialIndex;
currTrial.currentOffers	= currentOffers;
currTrial.chosenOffer	= [];
currTrial.gotjuice		= 0;
% currTrial.C			= C;
TrialRecord.User.currTrial	= currTrial;
TrialRecord.User.trialRecord(trialNumber) = currTrial;
%
TrialRecord.User.blocksDone			= blocksDone;
TrialRecord.User.halfBlocksDone		= halfBlocksDone;
TrialRecord.User.trialsDoneThisHalfBlock = trialsDoneThisHalfBlock;

%
%output trial number to DAQ
if testtrial
    % 	output_trialnumber(-1, DaqInfo);
    DAQ.eventmarker(-1);
else
    % 	output_trialnumber(trialNumber);
    DAQ.eventmarker(trialNumber);
end


function [] = output_trialnumber(trialNumber, varargin)

persistent DaqDIO digoutflag z databits strobebit sbval numdatabits

if trialNumber == -1, %set trial-start time
    DAQ = varargin{1};
    digoutflag = 0;
    if isfield(DAQ.BehavioralCodes, 'DIO'),
        digoutflag = 1;
        DaqDIO = DAQ.BehavioralCodes.DIO;
        databits = DAQ.BehavioralCodes.DataBits.Index;
        databits = cat(2, databits{:});
        numdatabits = length(databits);
        strobebit = DAQ.BehavioralCodes.StrobeBit.Index;
        z = zeros(1, numdatabits+1);
        putvalue(DaqDIO, z);
    end
    sbval = DAQ.StrobeBitEdge - 1; %falling edge -> 0 or rising edge -> 1
    return
end

% Output trial number on digital port
if digoutflag,
    trialnum_14bits = dec2binvec(trialNumber, 14);
    trialnum_lobyte = [trialnum_14bits(1:7), 1];
    trialnum_hibyte = [trialnum_14bits(8:14), 1];
    
    bvec = z;
    bvec([databits strobebit]) = [trialnum_hibyte ~sbval];
    putvalue(DaqDIO, bvec);
    bvec(strobebit) = sbval;
    putvalue(DaqDIO, bvec);
    pause(25/1000);		%wait 25ms before sending low bit
    %
    bvec([databits strobebit]) = [trialnum_lobyte ~sbval];
    putvalue(DaqDIO, bvec);
    bvec(strobebit) = sbval;
    putvalue(DaqDIO, bvec);
end


