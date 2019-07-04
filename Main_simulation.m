clear all; close all; clc;

%% Initialize
% data = load('Data_Human_new.txt');
N = 10^3;
T = 0.1;
t = 0:T:(N-1)*T;

HRd = zeros(1,N)+1;   RRd = zeros(1,N)+1;
HRe = zeros(1,N);   RRe = zeros(1,N);
HRed = zeros(1,N);  RRed = zeros(1,N);
KH = zeros(1,N);    SPD = zeros(1,N);
HR = zeros(1,N);    RR = zeros(1,N);
Toq = zeros(1,N);   Leveld = zeros(1,N);
Level = zeros(1,N);

global Point; Point = [0 0.1 3 -3 0];

%% Reference

for k=1:N,
%     if k< N/5, HRd(k) = 1.2; RRd(k) = 1.1;
%     elseif k>=N/5 && k<N/3, HRd(k) = 1.4; RRd(k) = 1.3;
%     elseif k>=N/3 && k<N/2, HRd(k) = 2; RRd(k) = 1.6;    
%     else, HRd(k) = 2.2; RRd(k) = 1.7;
%     end
    if k< N/4, HRd(k) = 1; RRd(k) = 1; end
    if k>=N/4 && k<N/2, HRd(k) = 2.0; RRd(k) = 1.6; end
    
    if k<N/3.8, Leveld(k) = 1; end
    if k>=N/3.8 && k<N/2, Leveld(k) = 5; end
    
    if k>=N/2 && k<N/1.5, HRd(k) = 1.4; RRd(k) = 1.3; Leveld(k) = 3; end
    if k>=N/1.5, HRd(k) = 1.8; RRd(k) = 1.45; Leveld(k) = 4;end
    
    
%     if k>=N/3 && k<N/2.5, Toq(k) = -0.3*(1-cos(2*pi*(k-N/3)/(N/2.5-N/3))); end
%     if k>=N/2 && k<N/1.8, Toq(k) = -0.2*(1-cos(2*pi*(k-N/2)/(N/1.8-N/2))); end
    
end

%% Fuzzifier Input
% Heart rate (beat per minute):
% Normally:  60 ~ 100 beat per minute
% Walking: 110 ~ 120
% Maximum: 200 - age
% HR = linspace(1,2.3,N);
% HR_mbs0 = linspace(1.2,2.1,5); % NH NL ZE PL PH
% HRd = linspace(1,2.3,N);
% HRd_mbs0 = linspace(1.2,2.1,5); % NH NL ZE PL PH
HRe_mbs = linspace(-0.1,0.1,N);
HRe_mbs0 = linspace(-0.05,0.05,5); % NH NL ZE PL PH
HRed_mbs = linspace(-3,3,N);
HRed_mbs0 = linspace(-1.5,1.5,3); % NH NL ZE PL PH

% Respiration rate (breaths per minute):
% 12 ~ 20 : Normal
% <12 : Low, > 25: High
% RR = linspace(1,1.8,N);
% RR_mbs0 = linspace(1.3,1.6,5);  % NH NL ZE PL PH
% RRd = linspace(1,1.8,N);
% RRd_mbs0 = linspace(1.3,1.6,5);  % NH NL ZE PL PH
RRe_mbs = linspace(-0.1,0.1,N);
RRe_mbs0 = linspace(-0.05,0.05,5);  % NH NL ZE PL PH
RRed_mbs = linspace(-3,3,N);
RRed_mbs0 = linspace(-1.5,1.5,3);  % NH NL ZE PL PH

% Knee Torque in Swing
Toq2 = linspace(-1,1,N);
Toq2_mbs0 = linspace(-0.75,0.75,5);  % NH NL ZE PL PH

%% Fuzzifier Output
KH_mbs  = linspace(-0.3,0.3,N);          % Knee Height: 1 ~ 5 level
KH_mbs0 = linspace(-0.15,0.15,5);          % NH, NL, ZN, ZP, PL, PH

SPD_mbs  = linspace(-0.4,0.4,N);     % Speed : 0.5 ~ 2.5 km/h
SPD_mbs0 = linspace(-0.2,0.2,5);     % NH, NL, ZN, ZP, PL, PH


HRlv_mbs = linspace(0,2.5,N);
HRlv_mbs0 = linspace(1,2,5); % NH NL ZE PL PH

RRlv_mbs = linspace(0,1.8,N);
RRlv_mbs0 = linspace(1,1.6,5); % NH NL ZE PL PH

Level_mbs  = linspace(0,6,N);     % Speed : 0.5 ~ 2.5 km/h
Level_mbs0 = linspace(1,5,5);     % NH, NL, ZN, ZP, PL, PH



%% Membership function
% Input
[HRe_mbsfn, HRe_sum, sHRe] = Gauss_mbs(HRe_mbs,HRe_mbs0,N);
[RRe_mbsfn, RRe_sum, sRRe] = Gauss_mbs(RRe_mbs,RRe_mbs0,N);
[HRed_mbsfn, HRed_sum, sHRed] = Gauss_mbs(HRed_mbs,HRed_mbs0,N);
[RRed_mbsfn, RRed_sum, sRRed] = Gauss_mbs(RRed_mbs,RRed_mbs0,N);
[Toq2_mbsfn, Toq2_sum, sToq2] = Gauss_mbs(Toq2,Toq2_mbs0,N);

% Output
[KH_mbsfn, KH_sum, sKH] = Gauss_mbs(KH_mbs,KH_mbs0,N);
[SPD_mbsfn, SPD_sum, sSPD] = Gauss_mbs(SPD_mbs,SPD_mbs0,N);

% Level
[HRlv_mbsfn, HRlv_sum, sHRlv] = Gauss_mbs(HRlv_mbs,HRlv_mbs0,N);
[RRlv_mbsfn, RRlv_sum, sRRlv] = Gauss_mbs(RRlv_mbs,RRlv_mbs0,N);
[Level_mbsfn, Level_sum, sLevel] = Gauss_mbs(Level_mbs,Level_mbs0,N);


while(1) % Discrete Check
    % Point
    [Point_HRe_mbsfn, Point_HRe_sum, sHRe] = Gauss_mbs(Point(1),HRe_mbs0,1);
    [Point_RRe_mbsfn, Point_RRe_sum, sRRe] = Gauss_mbs(Point(2),RRe_mbs0,1);
    [Point_HRed_mbsfn, Point_HRed_sum, sHRed] = Gauss_mbs(Point(3),HRed_mbs0,1);
    [Point_RRed_mbsfn, Point_RRed_sum, sRRed] = Gauss_mbs(Point(4),RRed_mbs0,1);
    [Point_Toq2_mbsfn, Point_Toq2_sum, sToq2] = Gauss_mbs(Point(5),Toq2_mbs0,1);
    
    [num,sum_num,mu,mu_num] = fuzzyrule_specific(Point_HRe_mbsfn,Point_RRe_mbsfn,Point_HRed_mbsfn,Point_RRed_mbsfn,Point_Toq2_mbsfn,KH_mbs0,SPD_mbs0);
    KH_num  = num(1);  KH_mu  = mu(1);
    SPD_num = num(2);  SPD_mu = mu(2);
    [Point_KH_mbsfn, Point_KH_sum, sKH] = Gauss_mbs(KH_mu,KH_mbs0,1);
    [Point_SPD_mbsfn, Point_SPD_sum, sSPD] = Gauss_mbs(SPD_mu,SPD_mbs0,1);
    fprintf('\nOutput : [KH = %1.4f (level), SPD = %1.4f (km/h)]\n',mu)
    %% [Membership function Input ] Figure
    figure('color','w')
    
    subplot(511);
    for k=1:length(HRe_mbs0), plot(HRe_mbs,HRe_mbsfn(:,k),'linewidth',2); hold on; end
    for k=1:length(HRe_mbs0), plot(Point(1),Point_HRe_mbsfn(:,k),'ro','linewidth',2); hold on; end
    %     title(sprintf('HR_{mbs}: [%1.4f  %1.4f  %1.4f  %1.4f  %1.4f]',Point_HR_mbsfn));
    ylabel('\mu_{ln1}(x)'); xlabel('bpm (beat / min)')
    
    subplot(512);
    for k=1:length(RRe_mbs0), plot(RRe_mbs,RRe_mbsfn(:,k),'linewidth',2); hold on; end
    for k=1:length(RRe_mbs0), plot(Point(2),Point_RRe_mbsfn(:,k),'ro','linewidth',2); hold on; end
    %     title(sprintf('RR_{mbs}: [%1.4f  %1.4f  %1.4f  %1.4f  %1.4f]',Point_RR_mbsfn));
    ylabel('\mu_{ln2}(x)'); xlabel('bpm (breaths / min)')
    
    subplot(513);
    for k=1:length(HRed_mbs0), plot(HRed_mbs,HRed_mbsfn(:,k),'linewidth',2); hold on; end
    for k=1:length(HRed_mbs0), plot(Point(3),Point_HRed_mbsfn(:,k),'ro','linewidth',2); hold on; end
    %     title(sprintf('HR_{mbs}: [%1.4f  %1.4f  %1.4f  %1.4f  %1.4f]',Point_HR_mbsfn));
    ylabel('\mu_{ln1}(x)'); xlabel('bpm rate (beat / min)')
    
    subplot(514);
    for k=1:length(RRed_mbs0), plot(RRed_mbs,RRed_mbsfn(:,k),'linewidth',2); hold on; end
    for k=1:length(RRed_mbs0), plot(Point(4),Point_RRed_mbsfn(:,k),'ro','linewidth',2); hold on; end
    %     title(sprintf('RR_{mbs}: [%1.4f  %1.4f  %1.4f  %1.4f  %1.4f]',Point_RR_mbsfn));
    ylabel('\mu_{ln2}(x)'); xlabel('bpm rate (breaths / min)')
    
    subplot(515);
    for k=1:length(Toq2_mbs0), plot(Toq2,Toq2_mbsfn(:,k),'linewidth',2); hold on; end
    for k=1:length(Toq2_mbs0), plot(Point(5),Point_Toq2_mbsfn(:,k),'ro','linewidth',2); hold on; end
    %     title(sprintf('Toq2_{mbs}: [%1.4f  %1.4f  %1.4f]',Point_Toq2_mbsfn));
    ylabel('\mu_{ln5}(x)'); xlabel('\tau_{nom} Nm/Nm')
    
    hFig = figure(1);
    set(hFig, 'Position', [300 100 500 800])
    set(gcf, 'renderer', 'painters');
    drawnow;
    %% [Membership function Input ] Figure
    figure('color','w')
    
    subplot(211);
    for k=1:length(KH_mbs0), plot(KH_mbs,KH_mbsfn(:,k),'linewidth',2); hold on; end
    for k=1:length(KH_mbs0), plot(KH_mu,Point_KH_mbsfn(:,k),'ro','linewidth',2); hold on; end
    %     title(sprintf('KH_{mbs}: [%1.4f  %1.4f  %1.4f]',Point_KH_mbsfn));
    ylabel('\mu_{out2}(x)'); xlabel('Knee height (Level)')
    
    subplot(212);
    for k=1:length(SPD_mbs0), plot(SPD_mbs,SPD_mbsfn(:,k),'linewidth',2); hold on; end
    for k=1:length(SPD_mbs0), plot(SPD_mu,Point_SPD_mbsfn(:,k),'ro','linewidth',2); hold on; end
    %     title(sprintf('SPD_{mbs}: [%1.4f  %1.4f  %1.4f %1.4f %1.4f]',Point_SPD_mbsfn));
    ylabel('\mu_{out3}(x)'); xlabel('Speed (km/h)')
    
    hFig = figure(2);
    set(hFig, 'Position', [800 400 500 450])
    set(gcf, 'renderer', 'painters');
    drawnow;
    break;
end

%% Fuzzy Control feedback loop

Level(1) = 1;
a = [-0.09453 0.4331 -0.04252 1];
b = [-0.09053 0.4031 -0.04052 1];

SPD(1) = 0.5; % 0.5 ~ 2.5
KH(1) = 1;    % 1 ~ 5
k=1;
HR(k) = (a(1)*(SPD(k)-0.5).^3 + a(2)*(SPD(k)-0.5).^2 + a(3)*(SPD(k)-0.5) + a(4)) .* (0.5*tanh(1.6*(KH(k)-2.5))+1.492) * 0.5 + 0.5;
RR(k) = (b(1)*(SPD(k)-0.5).^3 + b(2)*(SPD(k)-0.5).^2 + b(3)*(SPD(k)-0.5) + b(4)) .* ((0.25*tanh(1.1*(KH(k)-2.7)))+1.2384) * 0.5 + 0.5;
  
j=1;
for k=2:N
    
    HR(k) = (a(1)*(SPD(k-1)-0.5).^3 + a(2)*(SPD(k-1)-0.5).^2 + a(3)*(SPD(k-1)-0.5) + a(4)) .* (0.5*tanh(1.6*(KH(k-1)-2.5))+1.492) * 0.5 + 0.5;
    RR(k) = (b(1)*(SPD(k-1)-0.5).^3 + b(2)*(SPD(k-1)-0.5).^2 + b(3)*(SPD(k-1)-0.5) + b(4)) .* ((0.25*tanh(1.1*(KH(k-1)-2.7)))+1.2384) * 0.5 + 0.5;
    
    HRe(k) = HRd(k) - HR(k);
    RRe(k) = RRd(k) - RR(k);
    HRed(k) = (HRe(k)-HRe(k-1))/T;
    RRed(k) = (RRe(k)-RRe(k-1))/T;
    
    % Membership function
    [Point_HRe_mbsfn, Point_HRe_sum] = Gauss_mbs_sinput(HRe(k),HRe_mbs0,1,sHRe);
    [Point_RRe_mbsfn, Point_RRe_sum] = Gauss_mbs_sinput(RRe(k),RRe_mbs0,1,sRRe);
    [Point_HRed_mbsfn, Point_HRed_sum] = Gauss_mbs_sinput(HRed(k),HRed_mbs0,1,sHRed);
    [Point_RRed_mbsfn, Point_RRed_sum] = Gauss_mbs_sinput(RRed(k),RRed_mbs0,1,sRRed);
    [Point_Toq2_mbsfn, Point_Toq2_sum] = Gauss_mbs_sinput(Toq(k),Toq2_mbs0,1,sToq2);
    
    if k==j*2
    % Funzzy Rule
    [num,sum_num,mu,mu_num] = fuzzyrule_specific(Point_HRe_mbsfn,Point_RRe_mbsfn,Point_HRed_mbsfn,Point_RRed_mbsfn,Point_Toq2_mbsfn,KH_mbs0,SPD_mbs0);
    j=j+1;
    else, mu=[0 0];
    end
    KH(k)  = -mu(1) + KH(k-1);
    SPD(k) = -mu(2) + SPD(k-1);
    
    % Level
    [Point_HRlv_mbsfn, Point_HRlv_sum] = Gauss_mbs_sinput(HR(k),HRlv_mbs0,1,sHRlv);    
    [Point_RRlv_mbsfn, Point_RRlv_sum] = Gauss_mbs_sinput(RR(k),RRlv_mbs0,1,sRRlv);
    [numlevel,sum_numlevel,mulevel,mu_numlevel] = fuzzyrule_level(Point_HRlv_mbsfn,Point_RRlv_mbsfn,Level_mbs0);
    Level(k) = mulevel;
    
    % Saturation of control input
    if SPD(k) > 2.5, SPD(k) = 2.5; elseif SPD(k) < 0.5, SPD(k) = 0.5; end
    if KH(k) > 5, KH(k) = 5; elseif KH(k) < 1, KH(k) = 1; end
end

% % Figure;
figure('color','w');
subplot(511); 
plot(t,Leveld,'r','linewidth',2); hold on;
plot(t,Level,'b:','linewidth',2); hold on; ylabel('Level')
legend('Desired Level','Actual Level')
axis([0 100 0.5 5]);
set(gca,'fontsize',12,'YTick',[1 3 5],'box','off');
xlabel('Session (%)'); grid on;

subplot(512); 
plot(t,HR,'b','linewidth',2); ylabel('Amplitude');
xlabel('Session (%)');
set(gca,'fontsize',12,'YTick',[1 1.5 2],'box','off');
axis([0 100 1 2]);
grid on;

subplot(513); 
plot(t,RR,'b','linewidth',2); ylabel('Amplitude');
xlabel('Session (%)')
set(gca,'fontsize',12,'YTick',[1 1.5 2],'box','off');
axis([0 100 1 2]);
grid on;

subplot(514); plot(t,SPD,'b','linewidth',2); hold on; ylabel('£ö (km/h)');
xlabel('Session (%)')
set(gca,'fontsize',12,'YTick',0.5:1:2.5,'box','off');
axis([0 100 0.5 2.5]);
grid on;

subplot(515); plot(t,KH,'b','linewidth',2); hold on; ylabel('Level');
xlabel('Session (%)')
set(gca,'fontsize',12,'YTick',[1 2 3],'box','off');
axis([0 100 1 3]);
grid on;

hFig = figure(3);
set(hFig, 'Position', [300 100 1000 600])
set(gcf, 'renderer', 'painters');
drawnow;
