clear parsession
filename = ['TT_' , session];	
eval([filename])

monkey = parsession.monkey;
day = parsession.day;
sessnum = parsession.sessnum;
% 
dirroot = ['C:\Experiments\TwoTasks\Data\',parsession.mnkdir, filesep, day, filesep];

fileroot_ML = ['TT', session(2:end)];
fileroot_CED = session(1:end-1);
filenums_CED = parsession.CEDfilenums;

cells_td = parsession.clusters;


uprobe = parsession.uprobe;
elnum = parsession.elnum;
channels = parsession.locations(:,1);
channels_label = parsession.locations(:,end); % 4 for OFC; 9 for dLPFC; 10 for vLPFC
channels_location = parsession.locations(:,1:end-1);

%flags
flags.trialstart	= 21;
flags.prefixon      = 24; 
flags.fixon			= 25; % cue for trial type (JC or SO)
flags.fixoff		= 27;
flags.prefixoff     = 28; 
flags.offeron		= 30;
flags.offeroff		= 31;
flags.offer2on		= 32;
flags.offer2off		= 33;
flags.sacctgton		= 35;
flags.sacctgtoff	= 36; 
flags.choicemade	= 39;
flags.juice1		= 43;
flags.juice2		= 44;
flags.juice3		= 45;
flags.trialend		= 60;
flags.outcome		= 50; 

flags.poss_outcomes = [flags.juice1, flags.juice2, flags.juice3];


