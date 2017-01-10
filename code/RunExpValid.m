%% Experimental Validation of the Proposed Reference-free Delamination Detection
%
%% Contact
%  Name  : Chul Min Yeum
%  Email : chulminy@gmail.com
%  Please contact me if you have a question or find a bug. You can use a
%  "issues" page in the github.

%% Description
%
% This code is used for the following paper: Chul Min Yeum, Hoon Sohn,
% Hyung Jin Lim, and Jeong Beom Ihn. “Reference-Free Delamination Detection
% Using Lamb Waves.” Structural Control and Health Monitoring 21, no. 5
% (May 1, 2014): 675–84.
%
% This code is tested under only window machine so you may find errors in
% other OS, but it may be minor errors related with setting up file pathes
% or toolbox download in "Parameters". I am open to assist you to work this
% code on differt OS.
%
% You can feel free to modify this code but I would recommend cite my
% paper.

%% Reference
% 
% Chul Min Yeum, Hoon Sohn, Hyung Jin Lim, and Jeong Beom Ihn.
% “Reference-Free Delamination Detection Using Lamb Waves.” Structural
% Control and Health Monitoring 21, no. 5 (May 1, 2014): 675–84.
%
% Chul Min Yeum, Hoon Sohn, Jeong Beom Ihn, Hyung Jin Lim, “Reference-free
% delamination detection using Lamb wave time delay,” the 8h International
% Workshop on Structural Health Monitoring, Stanford, CA, September 13-15,
% 2011.

%% Setup
clear; clc; close all; format shortg; warning off;

% folder setup
folderBase  = fullfile(cd(cd('..')),'data'); 
folderExp   = fullfile(folderBase,'exp');
folderPrc   = fullfile(folderBase,'prc_data'); clearvars folderBase;
folderOut   = fullfile(cd(cd('..')),'output'); 

% gabor function
dt      = 2e-8;         % sampling frequency
fc      = 80*1000;      % driving frequency (80k)
sigma   = 7/(2*pi*fc);  % 7: number of cycle.

% All signals are cropped so that they have 'cropTp' time points, but they
% are different intial time points. The signal should not include
% refelection.
cropTp  = 13000;

% time range: we crop the signals having 'cropTp' time points
x_time  = (0:cropTp-1)*dt; 

% S0 range (to extract a A0 mode from each signal)
S0Range = 1:3000;

%% Decomposition of A0 signals
if ~exist(fullfile(folderPrc,'DataA0.mat'),'file');
    LoadSignal;
else
    load(fullfile(folderPrc,'DataA0.mat'));
end

%% Extract waveforms
% I can't publish a code for matching pursuit due to copyright. Here, I
% make a code with a similar concept.
if ~exist(fullfile(folderPrc,'DataM.mat'),'file');
    %% Intact and Temp 20
    M0          = A0IntT20;
    xdata       = 0:dt:(numel(M0)-1)*dt;
    
    [amp,u]  = max(abs(M0));
    x0  = [amp fc u*dt sigma 0];
    lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,   -pi];
    ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1,    pi];
    xM1 = lsqcurvefit(@funGaussModSine,x0,xdata,M0,lb, ub);
    M1  = funGaussModSine(xM1,xdata); clearvars x0 lb ub amp u
    
    tmp = M0-M1; 
    M2  = zeros(1,cropTp);
    xM2 = xM1; xM2(3) = xM2(3)-dt;
    while xM1(3)>xM2(3)
        [amp,u]  = max(abs(tmp-M2));
        x0  = [amp fc u*dt sigma 0];
        lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,-pi];
        ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1, pi];
        xM2 = lsqcurvefit(@funGaussModSine,x0,xdata,M0-M1,lb, ub);
        M2  = funGaussModSine(xM2,xdata); 
    end; clearvars xM2 x0 lb ub amp u
    
    residual   = M0-M1-M2;
    errorBound = std(residual)*3;
    
    DataInt20.M0 = M0;
    DataInt20.M1 = M1;
    DataInt20.M2 = M2;
    DataInt20.errorBound = errorBound;
    
    %% Intact and Temp 50
    M0          = A0IntT50;
    xdata       = 0:dt:(numel(M0)-1)*dt;
    
    [amp,u]  = max(abs(M0));
    x0  = [amp fc u*dt sigma 0];
    lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,   -pi];
    ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1,    pi];
    xM1 = lsqcurvefit(@funGaussModSine,x0,xdata,M0,lb, ub);
    M1  = funGaussModSine(xM1,xdata); clearvars x0 lb ub amp u
    
    tmp = M0-M1; 
    M2  = zeros(1,cropTp);
    xM2 = xM1; xM2(3) = xM2(3)-dt;
    while xM1(3)>xM2(3)
        [amp,u]  = max(abs(tmp-M2));
        x0  = [amp fc u*dt sigma 0];
        lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,-pi];
        ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1, pi];
        xM2 = lsqcurvefit(@funGaussModSine,x0,xdata,M0-M1,lb, ub);
        M2  = funGaussModSine(xM2,xdata); 
    end; clearvars xM2 x0 lb ub amp u
    
    residual   = M0-M1-M2;
    errorBound = std(residual)*3;
    
    DataInt50.M0 = M0;
    DataInt50.M1 = M1;
    DataInt50.M2 = M2;
    DataInt50.errorBound = errorBound;
    
    %% Damage and Temp 20
    M0          = A0DmgT20;
    xdata       = 0:dt:(numel(M0)-1)*dt;
    
    [amp,u]  = max(abs(M0));
    x0  = [amp fc u*dt sigma 0];
    lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,   -pi];
    ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1,    pi];
    xM1 = lsqcurvefit(@funGaussModSine,x0,xdata,M0,lb, ub);
    M1  = funGaussModSine(xM1,xdata); clearvars x0 lb ub amp u
    
    tmp = M0-M1; 
    M2  = zeros(1,cropTp);
    xM2 = xM1; xM2(3) = xM2(3)-dt;
    while xM1(3)>xM2(3)
        [amp,u]  = max(abs(tmp-M2));
        x0  = [amp fc u*dt sigma 0];
        lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,-pi];
        ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1, pi];
        xM2 = lsqcurvefit(@funGaussModSine,x0,xdata,M0-M1,lb, ub);
        M2  = funGaussModSine(xM2,xdata); 
    end; clearvars xM2 x0 lb ub amp u
    
    residual   = M0-M1-M2;
    errorBound = std(residual)*3;
    
    DataDmg20.M0 = M0;
    DataDmg20.M1 = M1;
    DataDmg20.M2 = M2;
    DataDmg20.errorBound = errorBound;
    
    %% Damage and Temp 50
    M0          = A0DmgT50;
    xdata       = 0:dt:(numel(M0)-1)*dt;
    
    [amp,u]  = max(abs(M0));
    x0  = [amp fc u*dt sigma 0];
    lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,   -pi];
    ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1,    pi];
    xM1 = lsqcurvefit(@funGaussModSine,x0,xdata,M0,lb, ub);
    M1  = funGaussModSine(xM1,xdata); clearvars x0 lb ub amp u
    
    tmp = M0-M1; 
    M2  = zeros(1,cropTp);
    xM2 = xM1; xM2(3) = xM2(3)-dt;
    while xM1(3)>xM2(3)
        [amp,u]  = max(abs(tmp-M2));
        x0  = [amp fc u*dt sigma 0];
        lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,-pi];
        ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1, pi];
        xM2 = lsqcurvefit(@funGaussModSine,x0,xdata,M0-M1,lb, ub);
        M2  = funGaussModSine(xM2,xdata); 
    end; clearvars xM2 x0 lb ub amp u
    
    residual    = M0-M1-M2;
    errorBound  = std(residual)*3;
    
    DataDmg50.M0 = M0;
    DataDmg50.M1 = M1;
    DataDmg50.M2 = M2;
    DataDmg50.errorBound = errorBound;
    
    save(fullfile(folderPrc,'DataM.mat'),'DataInt20','DataInt50', ...
        'DataDmg20','DataDmg50');
else
    load(fullfile(folderPrc,'DataM.mat'),'DataInt20','DataInt50', ...
        'DataDmg20','DataDmg50');
end

%% Ploting outcome

% Load raw signal
if ~exist(fullfile(folderPrc,'DataRB.mat'),'file') && ...
        ~exist(fullfile(folderPrc,'DataDB.mat'),'file');
    LoadSignal;
else
    load(fullfile(folderPrc,'DataRB.mat'));
    load(fullfile(folderPrc,'DataDB.mat'));
end

%--------------------------------------------------------------------------
% Figure 1a
sigIdx     = 10; % among 2-;
h = figure('name','Figure1a'); 
plot(x_time*1000,RBIntT20(sigIdx,:)*1000,'b', ...
     x_time*1000,RBDmgT20(sigIdx,:)*1000,':r','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Undamaged (20^oC)','\bf Damaged (20^oC)'); legend('boxoff');
set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0','0.05','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-10 -5 0 5 10])
set(gca,'YTickLabel',{'-10','-5','0','5','10'})
xlabel('\bf Time(ms)');
ylabel('\bf Voltage(mV)');
ylim([-11 11]); xlim([0 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure1a'),'-djpeg','-r0')

% Figure 1b
h = figure('name','Figure1b');
normScaleY1 = max(DataInt20.M0);
normScaleY2 = max(DataDmg20.M0);
plot(x_time*1000,DataInt20.M0/normScaleY1,'b', ...
     x_time*1000,DataDmg20.M0/normScaleY2,':r','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Undamaged (20^oC)','\bf Damaged (20^oC)'); legend('boxoff');
set(gca,'XTick',[0.06 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0.06','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-1 -0.5 0 0.5 1])
set(gca,'YTickLabel',{'-1','-0.5','0','0.5','1'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitdue');
ylim([-1.2 1.2]); xlim([0.06 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure1b'),'-djpeg','-r0')

% Figure 1c
sigIdx     = 10; % among 2-;
h = figure('name','Figure1c'); 
plot(x_time*1000,RBIntT20(sigIdx,:)*1000,'b', ...
     x_time*1000,RBIntT50(sigIdx,:)*1000,':r','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Undamaged (20^oC)','\bf Undamaged (50^oC)'); legend('boxoff');
set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0','0.05','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-10 -5 0 5 10])
set(gca,'YTickLabel',{'-10','-5','0','5','10'})
xlabel('\bf Time(ms)');
ylabel('\bf Voltage(mV)');
ylim([-11 11]); xlim([0 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure1c'),'-djpeg','-r0')

% Figure 1d
h = figure('name','Figure1d');
normScaleY1 = max(DataInt20.M0);
normScaleY2 = max(DataInt50.M0);
plot(x_time*1000,DataInt20.M0/normScaleY1,'b', ...
     x_time*1000,DataInt50.M0/normScaleY2,':r','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Undamaged (20^oC)','\bf Undamaged (50^oC)'); legend('boxoff');
set(gca,'XTick',[0.06 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0.06','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-1 -0.5 0 0.5 1])
set(gca,'YTickLabel',{'-1','-0.5','0','0.5','1'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitdue');
ylim([-1.2 1.2]); xlim([0.06 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure1d'),'-djpeg','-r0')

%--------------------------------------------------------------------------

% Figure 3(a)
sigIdx     = 10; % among 2-;
h = figure('name','Figure3(a)'); 
plot(x_time*1000,RBDmgT20(sigIdx,:)*1000,'b', ...
     x_time*1000,DBDmgT20(sigIdx,:)*1000,':r','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Signal Ab','\bf Signal ab'); legend('boxoff');
set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0','0.05','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-10 -5 0 5 10])
set(gca,'YTickLabel',{'-10','-5','0','5','10'})
xlabel('\bf Time(ms)');
ylabel('\bf Voltage(mV)');
ylim([-11 11]); xlim([0 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure3(b)'),'-djpeg','-r0')

% Figure 3(b)
h = figure('name','Figure3(b)');
Sc = max(RBDmgT20(sigIdx,1:3000))/max(DBDmgT20(sigIdx,1:3000));
plot(x_time*1000,RBDmgT20(sigIdx,:)*1000,'b', ...
     x_time*1000,DBDmgT20(sigIdx,:)*1000*Sc,':r','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Signal Ab','\bf Signal ab x Sc'); legend('boxoff');
set(gca,'XTick',[0 0.05 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0','0.05','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-10 -5 0 5 10])
set(gca,'YTickLabel',{'-10','-5','0','5','10'})
xlabel('\bf Time(ms)');
ylabel('\bf Voltage(mV)');
ylim([-11 11]); xlim([0 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure3(b)'),'-djpeg','-r0')

% Figure 3(c)
h = figure('name','Figure3(c)');
sig = RBDmgT20(sigIdx,:)*1000 - DBDmgT20(sigIdx,:)*1000*Sc;
plot(x_time*1000,sig./max(sig),'b','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Signal Ab - Signal ab x Sc'); legend('boxoff');
set(gca,'XTick',[0.06 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0.06','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-1 -0.5 0 0.5 1])
set(gca,'YTickLabel',{'-1','-0.5','0','0.5','1'})
xlabel('\bf Time(ms)');
ylabel('\bf Voltage(mV)');
ylim([-1.2 1.2]); xlim([0.06 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure3(c)'),'-djpeg','-r0')

% Figure 4
h = figure('name','Figure4');
normScaleY = max(DataDmg20.M0);
plot(x_time*1000,DataDmg20.M0/normScaleY,'b', ...
     x_time*1000,DataDmg20.M1/normScaleY,':r', ...
     x_time*1000,DataDmg20.M2/normScaleY,'-.k','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Original','\bf A_0_,_T mode','\bf A_0_,_R mode'); 
legend('boxoff');
set(gca,'XTick',[0.06 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0.06','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-1 -0.5 0 0.5 1])
set(gca,'YTickLabel',{'-1','-0.5','0','0.5','1'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitdue');
ylim([-1.2 1.2]); xlim([0.06 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure4'),'-djpeg','-r0')

% Figure 9
h = figure('name','Figure9(a)');
plot(x_time*1000,DataInt20.M2,'b', ...
     x_time*1000,DataInt50.M2,':r','linewidth',3);
line([x_time(1)*1000 x_time(end)*1000], ...
    [DataInt20.errorBound DataInt20.errorBound],...
    'Color','b','linewidth',3); 
line([x_time(1)*1000 x_time(end)*1000], ...
    -[DataInt20.errorBound DataInt20.errorBound],...
    'Color','b','linewidth',3);
line([x_time(1)*1000 x_time(end)*1000], ...
  [DataInt50.errorBound DataInt50.errorBound],...
  'Color','r','lineStyle',':','linewidth',3);
line([x_time(1)*1000 x_time(end)*1000], ...
 -[DataInt50.errorBound DataInt50.errorBound],...
 'Color','r','lineStyle',':','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf M_2(20^oC)','\bf M_2(50^oC)'); legend('boxoff');
set(gca,'XTick',[0.06 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0.06','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-0.1 0 0.1])
set(gca,'YTickLabel',{'-0.1','0','0.1'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitdue');
ylim([-0.15 0.15]); xlim([0.06 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure9a'),'-djpeg','-r0')

h = figure('name','Figure9(b)');
plot(x_time*1000,DataDmg20.M2,'b', ...
     x_time*1000,DataDmg50.M2,':r','linewidth',3);
line([x_time(1)*1000 x_time(end)*1000], ...
    [DataDmg20.errorBound DataDmg20.errorBound],...
    'Color','b','linewidth',3); 
line([x_time(1)*1000 x_time(end)*1000], ...
    -[DataDmg20.errorBound DataDmg20.errorBound],...
    'Color','b','linewidth',3); 
line([x_time(1)*1000 x_time(end)*1000], ...
  [DataDmg50.errorBound DataDmg50.errorBound],...
  'Color','r','lineStyle',':','linewidth',3);
line([x_time(1)*1000 x_time(end)*1000], ...
 -[DataDmg50.errorBound DataDmg50.errorBound],...
 'Color','r','lineStyle',':','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf M_2(20^oC)','\bf M_2(50^oC)'); legend('boxoff');
set(gca,'XTick',[0.06 0.1 0.15 0.2 0.25])
set(gca,'XTickLabel',{'0.06','0.1','0.15','0.2','0.25'})
set(gca,'YTick',[-0.12 0 0.12])
set(gca,'YTickLabel',{'-0.12','0','0.12'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitdue');
ylim([-0.18 0.18]); xlim([0.06 0.25]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure9b'),'-djpeg','-r0')