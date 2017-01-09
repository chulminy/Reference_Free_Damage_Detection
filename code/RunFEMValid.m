%% FEM Validation of the Proposed Reference-free Delamination Detection
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
% other OS. I am open to assist you to work this code on differt OS.
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
folderFEM   = fullfile(folderBase,'fem');
folderPrc   = fullfile(folderBase,'prc_data'); clearvars folderBase;
folderOut   = fullfile(cd(cd('..')),'output'); 

% load simulation data
load(fullfile(folderFEM,'DataSim.mat'));

% Here is the description of data obtained from FEM simulation
%
% DBInt: (Exicitation) Disk at PZT A -> (Measurement) PZT B on intact plate
% RBInt: (Exicitation) Ring at PZT A -> (Measurement) PZT B on intact plate
% DBDmg: (Exicitation) Disk at PZT A -> (Measurement) PZT B on damage plate
% RBDmg: (Exicitation) Ring at PZT A -> (Measurement) PZT B on damage plate

% gabor function
dt      = 5e-7;         % sampling frequency
fc      = 100*1000;     % frequency (100k)
sigma   = 7/(2*pi*fc);  % 7: number of cycle.

% time range
x_time  = (0:350)*dt- 0.00004; % 0.00004 is triggering point.

%% Ploting outcome
% Comparison of the Lamb wave signals obtained from undamaged and damaged
% conditions (numerical simulation): The A0,R mode and the subsequent
% reflection are clearly observed if delamination exists.
% normalize the amplitude with maximum value from Ab

% Figure 6(a)
normScaleY = max(DBInt); 
h = figure('name','Figure6a'); plot(x_time*1000,DBInt/normScaleY,'b', ...
     x_time*1000,DBDmg/normScaleY,':r','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Undamaged','\bf Damage'); legend('boxoff');
set(gca,'XTick',[-0.03 0 0.03 0.06 0.09 0.12])
set(gca,'XTickLabel',{'-0.03','0','0.03','0.06','0.09','0.12'})
set(gca,'YTick',[-2 -1 0 1 2])
set(gca,'YTickLabel',{'-2','-1','0','1','2'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitude');
ylim([-1.2 1.2]); xlim([-0.03 0.135]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure6a'),'-djpeg','-r0')

% Remove S0 mode: 80 is determined based on the signal having only S0 mode.
A0Int = (DBInt./max(DBInt(1:80)))-(RBInt./max(RBInt(1:80)));
A0Dmg = (DBDmg./max(DBDmg(1:80)))-(RBDmg./max(RBDmg(1:80)));

% Figure 6(b)
normScaleY = max(A0Int); 
h = figure('name','Figure6b'); plot(x_time*1000,A0Int/normScaleY,'b', ...
     x_time*1000,A0Dmg/normScaleY,':r','linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Undamaged','\bf Damage'); legend('boxoff');
set(gca,'XTick',[-0.03 0 0.03 0.06 0.09 0.12])
set(gca,'XTickLabel',{'-0.03','0','0.03','0.06','0.09','0.12'})
set(gca,'YTick',[-2 -1 0 1 2])
set(gca,'YTickLabel',{'-2','-1','0','1','2'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitude');
ylim([-1.2 1.2]); xlim([-0.03 0.135]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure6b'),'-djpeg','-r0')

% Mode decomposition
M0          = A0Dmg/normScaleY;
xdata       = 0:dt:(numel(M0)-1)*dt;

[amp,u]  = max(abs(M0));
x0  = [amp fc u*dt sigma 0];
lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,   -pi];
ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1,    pi];
xM1 = lsqcurvefit(@funGaussModSine,x0,xdata,M0,lb, ub);
M1  = funGaussModSine(xM1,xdata); clearvars x0 lb ub amp u

[amp,u]  = max(abs(M0-M1));
x0  = [amp fc u*dt sigma 0];
lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,   -pi];
ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1,    pi];
xM2 = lsqcurvefit(@funGaussModSine,x0,xdata,M0-M1,lb, ub);
M2  = funGaussModSine(xM2,xdata); clearvars x0 lb ub amp u

[amp,u]  = max(abs(M0-M1-M2));
x0  = [amp fc u*dt sigma 0];
lb  = [ amp-0.1*amp,  fc-fc*0.1, (u-1000)*dt, sigma-sigma*0.1,   -pi];
ub  = [ amp+0.1*amp,  fc+fc*0.1, (u+1000)*dt, sigma+sigma*0.1,    pi];
xM3 = lsqcurvefit(@funGaussModSine,x0,xdata,M0-M1-M2,lb, ub);
M3  = funGaussModSine(xM3,xdata); clearvars x0 lb ub amp u

% Figure 7
normScaleY = max(M0); 
h = figure('name','Figure6b'); 
plot(x_time*1000,M1/normScaleY,'b', ...
     x_time*1000,M2/normScaleY,':r', ...
     x_time*1000,M3/normScaleY,'-.k', 'linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf Transmitted','\bf 1^s^t reflected','\bf 2^n^d reflected'); 
legend('location','northwest');legend('boxoff'); 
set(gca,'XTick',[-0.03 0 0.03 0.06 0.09 0.12])
set(gca,'XTickLabel',{'-0.03','0','0.03','0.06','0.09','0.12'})
set(gca,'YTick',[-2 -1 0 1 2])
set(gca,'YTickLabel',{'-2','-1','0','1','2'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitude');
ylim([-1.2 1.2]); xlim([-0.03 0.135]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Figure7'),'-djpeg','-r0')

% Extra
% plot four signals: DBInt,RBInt,DBDmg,RBDmg
normScaleY = max(RBInt); 
h = figure('name','Extra(Undamage)'); 
plot(x_time*1000,DBInt/normScaleY,'b', ...
     x_time*1000,RBInt/normScaleY,':r', 'linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf A(ring) to B','\bf A(disk) to B'); 
legend('boxoff'); 
set(gca,'XTick',[-0.03 0 0.03 0.06 0.09 0.12])
set(gca,'XTickLabel',{'-0.03','0','0.03','0.06','0.09','0.12'})
set(gca,'YTick',[-2 -1 0 1 2])
set(gca,'YTickLabel',{'-2','-1','0','1','2'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitude');
ylim([-1.2 1.2]); xlim([-0.03 0.135]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Extra(Undamage)'),'-djpeg','-r0')

normScaleY = max(RBDmg); 
h = figure('name','Extra(Damage)'); 
plot(x_time*1000,DBDmg/normScaleY,'b', ...
     x_time*1000,RBDmg/normScaleY,':r', 'linewidth',3);
set(gca,'fontsize',20,'linewidth',3,'fontweight','bold')
legend('\bf A(ring) to B','\bf A(disk) to B'); 
legend('boxoff'); 
set(gca,'XTick',[-0.03 0 0.03 0.06 0.09 0.12])
set(gca,'XTickLabel',{'-0.03','0','0.03','0.06','0.09','0.12'})
set(gca,'YTick',[-2 -1 0 1 2])
set(gca,'YTickLabel',{'-2','-1','0','1','2'})
xlabel('\bf Time(ms)');
ylabel('\bf Normalized amplitude');
ylim([-1.2 1.2]); xlim([-0.03 0.135]); 
set(h,'PaperPositionMode','auto','pos',[50 50 900 450])
print(fullfile(folderOut,'Extra(Damage)'),'-djpeg','-r0')

