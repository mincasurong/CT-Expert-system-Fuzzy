clear all; close all; clc;

% % European Journal of applied physiology (2006)
% % Ref : Physiological responses to nordic walking, walking and jogging
% x0 = [0 0.2 0.3 0.4 0.5 0.6 1.2 1.3 1.5 1.65 1.8 2.1 2.4 2.5 2.7 2.8 2.9 3 3.5 4]*3.6; % km/h, walking speed for normal people
% x = x0 / 4;            % For patients
% HR = [75 75 76 76 78 80 95 100 107 115 125 150 152 155 156 156.5 157 157 157 158]/75; % bpm, beats per minute
% 
% figure; plot(x,HR)
% 
% %% Fit: 'untitled fit 1'.
% [xData, yData] = prepareCurveData( x, HR );
% 
% % Set up fittype and options.
% ft = fittype( 'poly4' );
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'HR vs. x', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel x
% ylabel HR
% grid on
% 

N=50; 
x=linspace(0.5,2.5,N);
%% HR vs spd
a = [-0.09453 0.4331 -0.04252 1];
for k=1:N;
    HR1(k) = 0.5*(a(1)*(x(k)-0.5)^3 + a(2)*(x(k)-0.5)^2 + a(3)*(x(k)-0.5)) + a(4);
end

%% RR vs spd
% for k=1:N;
%     if x(k)<0.5, RR1(k) = 0.1/0.5*x(k)+1; end
%     if x(k)>=0.5 && x(k)<2.5, RR1(k) = 1/2*(x(k)-0.5) + 1.1; end
%     if x(k)>=2.5, RR1(k) = 0.1/0.5*(x(k)-2.5) + 2.1; end
% end
a = [-0.09053 0.4031 -0.04052 1];
for k=1:N;
    RR1(k) = 0.5*(a(1)*(x(k)-0.5)^3 + a(2)*(x(k)-0.5)^2 + a(3)*(x(k)-0.5)) + a(4);
end

%% Figure
figure('color','w');
subplot(211); plot(x,HR1,'b'); ylabel('HR (beats/min)') 
subplot(212); plot(x,RR1,'r'); ylabel('RR (breaths/min)')
xlabel('Speed (km/h)');

hFig = figure(1);
set(hFig, 'Position', [300 100 400 400])
set(gcf, 'renderer', 'painters');
drawnow;
%% HR vs KH
y=linspace(1,5,N);
HR2=(0.25*tanh(1.6*(y-2.5)))+1.2459;

%% RR vs KH
RR2=(0.25*tanh(1.1*(y-2.7)))+1.2384;

figure('color','w');
subplot(211); plot(y,HR2,'b'); ylabel('HR (beats/min)') 
subplot(212); plot(y,RR2,'r'); ylabel('RR (breaths/min)')
xlabel('KH (level)');

hFig = figure(2);
set(hFig, 'Position', [300 100 400 400])
set(gcf, 'renderer', 'painters');
drawnow;

%% Bio vs spd vs KH
a = [-0.09453 0.4331 -0.04252 1];
b = [-0.09053 0.4031 -0.04052 1];
for k=1:N;
    HR(k) = a(1)*x(k)^3 + a(2)*x(k)^2 + a(3)*x(k) + a(4) + 0.5*tanh(1.6*(y(k)-2.5))+1.492;
    if x(k)<0.5, RR(k) = 0.1/0.5*x(k)+1; end
    if x(k)>=0.5 && x(k)<2.5, RR(k) = 1/2*(x(k)-0.5) + 1.1; end
    if x(k)>=2.5, RR(k) = 0.1/0.5*(x(k)-2.5) + 2.1; end
    RR(k) = RR(k) + 0.5*tanh(1.1*(y(k)-2.7))+1.477;
end

[X,Y] = meshgrid(x,y);
HR = (a(1)*(X-0.5).^3 + a(2)*(X-0.5).^2 + a(3)*(X-0.5) + a(4)) .* (0.5*tanh(1.6*(Y-2.5))+1.492) * 0.5 + 0.5;
RR = (b(1)*(X-0.5).^3 + b(2)*(X-0.5).^2 + b(3)*(X-0.5) + b(4)) .* ((0.25*tanh(1.1*(Y-2.7)))+1.2384) * 0.5 + 0.5;
figure('color','w');
subplot(211); s=surf(X,Y,HR); s.EdgeColor = 'none';  colorbar
axis([0.5 2.5 1 5 1 2.5])
xlabel('Speed (km/h)'); ylabel('KH (level)'); zlabel('HR (beats/min)') 
subplot(212); s=surf(X,Y,RR); s.EdgeColor = 'none'; colorbar
axis([0.5 2.5 1 5 1 2])
xlabel('Speed (km/h)'); ylabel('KH (level)'); zlabel('RR (breaths/min)')
maxHR = 2.383; 
maxRR = 1.842;

%% Plant
HRtest = @(q,w) (a(1)*(q-0.5).^3 + a(2)*(q-0.5).^2 + a(3)*(q-0.5) + a(4)) .* (0.5*tanh(1.6*(w-2.5))+1.492) * 0.5 + 0.5;
RRtest = @(q,w) (b(1)*(q-0.5).^3 + b(2)*(q-0.5).^2 + b(3)*(q-0.5) + b(4)) .* ((0.25*tanh(1.1*(w-2.7)))+1.2384) * 0.5 + 0.5;

q = input('Speed (0.5 ~ 2.5 km/h) = '); w = input('Knee height (1~5 level) = ');
Answer = [HRtest(q,w), RRtest(q,w)];
fprintf('\n HR = %1.4f, RR = %1.4f \n',Answer)