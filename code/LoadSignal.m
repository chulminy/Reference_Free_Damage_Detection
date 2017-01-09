%% Experimental setup

% Reference 

% [1] Chul Min Yeum, Hoon Sohn, Hyung Jin Lim, and Jeong Beom
% Ihn. “Reference-Free Delamination Detection Using Lamb Waves.” Structural
% Control and Health Monitoring 21, no. 5 (May 1, 2014): 675–84.
%
% [2] Chul Min Yeum, Hoon Sohn, Jeong Beom Ihn, Hyung Jin Lim,
% “Reference-free delamination detection using Lamb wave time delay,” the
% 8h International Workshop on Structural Health Monitoring, Stanford, CA,
% September 13-15, 2011.

% [1]: Journal paper, [2]: Conference paper

% The original PZT configuraition is below. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Dual PZT A  %%%%%%%%%%%%%%  PZT B  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please take a look at the folder "data-exp"
% Temp_20 and Temp_50 are signals measured under 0 celcius, room
% temperature, and 50 celcius, respectively.

% Each set of signals, triggering time points are different. 

%% Data: Intact and Temp 20
folder = fullfile(folderExp,'intact','Temp_20');
A0 = zeros(20,cropTp);
RB = zeros(20,cropTp);DB = zeros(20,cropTp);
for jj=1:20
    load(fullfile(folder,['80k_' int2str(jj), 'k_RB'])); 
    eval(['RBTmp = X80k_' int2str(jj), 'k_RB(10000:10000+cropTp-1);']); 
    load(fullfile(folder,['80k_' int2str(jj), 'k_DB'])); 
    eval(['DBTmp = X80k_' int2str(jj), 'k_DB(10000:10000+cropTp-1);']); 
    A0(jj,:) = (RBTmp./max(RBTmp(S0Range)) - ...
                DBTmp./max(DBTmp(S0Range)))';clearvars X*;
    RB(jj,:) = RBTmp;
    DB(jj,:) = DBTmp;
end; clearvars RBTmp DBTmp

% Peaks of 20 signals are slightly unmatched. Following process is to match
% a peck of each signal.
peakLoc = zeros(20,1);
for jj=1:20
    [~,I] = max(A0(jj,:)); 
    peakLoc(jj) = I;
end
peakDiff = ceil(median(peakLoc)-peakLoc);
for jj=1:20
    A0(jj,:) = circshift(A0(jj,:)', peakDiff(jj))';
end; clearvars peakDiff peakLoc;

A0IntT20 = mean(A0); clearvars A0;
RBIntT20 = RB; clearvars RB;
DBIntT20 = DB; clearvars DB;

%% Data: Intact and Temp 50
folder = fullfile(folderExp,'intact','Temp_50');
A0 = zeros(20,13000);
for jj=1:20
    load(fullfile(folder,['80k_' int2str(jj), 'k_RB'])); 
    eval(['RBTmp = X80k_' int2str(jj), 'k_RB(25000:25000+cropTp-1);']); 
    load(fullfile(folder,['80k_' int2str(jj), 'k_DB'])); 
    eval(['DBTmp = X80k_' int2str(jj), 'k_DB(25000:25000+cropTp-1);']); 
    A0(jj,:) = (RBTmp./max(RBTmp(S0Range)) - ...
                DBTmp./max(DBTmp(S0Range)))';clearvars X*;
    RB(jj,:) = RBTmp;
    DB(jj,:) = DBTmp;
end; clearvars RBTmp DBTmp

% Peaks of 20 signals are slightly unmatched. Following process is to match
% a peck of each signal.
peakLoc = zeros(20,1);
for jj=1:20
    [~,I] = max(A0(jj,:)); 
    peakLoc(jj) = I;
end
peakDiff = ceil(median(peakLoc)-peakLoc);
for jj=1:20
    A0(jj,:) = circshift(A0(jj,:)', peakDiff(jj))';
end; clearvars peakDiff peakLoc;

A0IntT50 = mean(A0); clearvars A0;
RBIntT50 = RB; clearvars RB;
DBIntT50 = DB; clearvars DB;

%% Data: Damage and Temp 20
folder = fullfile(folderExp,'damage','Temp_20');
A0 = zeros(20,13000);
for jj=1:20
    load(fullfile(folder,['80k_' int2str(jj), 'k_RB'])); 
    eval(['RBTmp = X80k_' int2str(jj), 'k_RB(50000:50000+cropTp-1);']); 
    load(fullfile(folder,['80k_' int2str(jj), 'k_DB'])); 
    eval(['DBTmp = X80k_' int2str(jj), 'k_DB(50000:50000+cropTp-1);']); 
    A0(jj,:) = (RBTmp./max(RBTmp(S0Range)) - ...
                DBTmp./max(DBTmp(S0Range)))';clearvars X*;
    RB(jj,:) = RBTmp;
    DB(jj,:) = DBTmp;
end; clearvars RBTmp DBTmp

% Peaks of 20 signals are slightly unmatched. Following process is to match
% a peck of each signal.
peakLoc = zeros(20,1);
for jj=1:20
    [~,I] = max(A0(jj,:)); 
    peakLoc(jj) = I;
end
peakDiff = ceil(median(peakLoc)-peakLoc);
for jj=1:20
    A0(jj,:) = circshift(A0(jj,:)', peakDiff(jj))';
end; clearvars peakDiff peakLoc;

A0DmgT20 = mean(A0); clearvars A0;
RBDmgT20 = RB; clearvars RB;
DBDmgT20 = DB; clearvars DB;

%% Data: Damage and Temp 50
folder = fullfile(folderExp,'damage','Temp_50');
A0 = zeros(20,13000);
for jj=1:20
    load(fullfile(folder,['80k_' int2str(jj), 'k_RB'])); 
    eval(['RBTmp = X80k_' int2str(jj), 'k_RB(50000:50000+cropTp-1);']); 
    load(fullfile(folder,['80k_' int2str(jj), 'k_DB'])); 
    eval(['DBTmp = X80k_' int2str(jj), 'k_DB(50000:50000+cropTp-1);']); 
    A0(jj,:) = (RBTmp./max(RBTmp(S0Range)) - ...
                DBTmp./max(DBTmp(S0Range)))';clearvars X*;
    RB(jj,:) = RBTmp;
    DB(jj,:) = DBTmp;
end; clearvars RBTmp DBTmp

% Peaks of 20 signals are slightly unmatched. Following process is to match
% a peck of each signal.
peakLoc = zeros(20,1);
for jj=1:20
    [~,I] = max(A0(jj,:)); 
    peakLoc(jj) = I;
end
peakDiff = ceil(median(peakLoc)-peakLoc);
for jj=1:20
    A0(jj,:) = circshift(A0(jj,:)', peakDiff(jj))';
end; clearvars peakDiff peakLoc;

A0DmgT50 = mean(A0); clearvars A0;
RBDmgT50 = RB; clearvars RB;
DBDmgT50 = DB; clearvars DB;

save(fullfile(folderPrc,'DataA0.mat'),'A0IntT20','A0IntT50', ...
    'A0DmgT20','A0DmgT50')
save(fullfile(folderPrc,'DataRB.mat'),'RBIntT20','RBIntT50', ...
    'RBDmgT20','RBDmgT50')
save(fullfile(folderPrc,'DataDB.mat'),'DBIntT20','DBIntT50', ...
    'DBDmgT20','DBDmgT50')
