%% Main

clear all
close all
clc

Npop=75;


%% Episodes list

% 1st temporal sequence)
ep1 = [1 2 3 4];
ep2 = [71 70 7 8 9 10];
ep3 = [11 12 13 14 15];
ep4 = [16 17 18 19 20 21];
ep5 = [22 23 24 25];

% 2nd temporal sequence)
ep6 = [26 27 28 29];
ep7 = [30 13 32 33 34 35];
ep8 = [36 37 38 39 40];
ep9 = [41 42 43 44 45 46];
ep10 = [47 48 49 50];

% 3rd temporal sequence)
ep11 = [51 52 53 54];
ep12 = [55 13 57 58 59 60];
ep13 = [61 62 44 64 65];
ep14 = [66 67 68 69 70 71];
ep15 = [72 73 74 75];


%% Simulation parameters

t_sim=1.45;   %simulation length
dt=0.0001;    %integration step
t=0:dt:t_sim;
T=length(t);

%load trained synapses
% load <<INSERT YOUR DIRECTORY HERE>> ... synapses_SHARED

%fixed synapses between layers
Wp_L1WM=eye(Npop)*100; % WM-->L1
Wp_WML1=eye(Npop)*0;   % L1-->WM
Wp_L2L1=eye(Npop)*300; % L1-->L2

%delays between layers
D_intraWM=1;
D_intraL1=1;
D_WML1=300;
D_L1WM=D_WML1;
D_thetaL1=100;
D_L1L2=10;
D_L2L1=D_L1L2;


%% Simulation

INPUT_WM=zeros(Npop,length(t));
buff=zeros(Npop, round(0.05/dt));
pos1=[1]; %first features of episode 1
buff(pos1,:)=1;
INPUT_WM(:,51:550)=buff;
buff=zeros(Npop, round(0.05/dt));
pos2=[26]; %first features of episode 6 (1st ep of the 2nd sequence)
buff(pos2,:)=1;
INPUT_WM(:,4301:4800)=buff;
buff=zeros(Npop, round(0.05/dt));
pos3=[51]; %first features of episode 11 (1st ep of the 3rd sequence)
buff(pos3,:)=1;
INPUT_WM(:,8701:9200)=buff;

RETRIEVAL_sim


%%  Figures

figure (1)

subplot(311), title('Working Memory'), hold on, ylabel('Zp (Hz)'), xlabel('time (s)')
plot(t,(zp0(1,:)),'b','LineWidth',0.75)
plot(t,(zp0(26,:)),'b--','LineWidth',0.75)
plot(t,(zp0(51,:)),'b:','LineWidth',0.75)
plot(t,(zp0([2:25 27:50 52:end],:)),'b:','LineWidth',0.75)
legend('Input ep1','Input ep2','Input ep3','Position',[0.913020854111298,0.741101223163884,0.07265625,0.07174638487208])
xlim([0 length(t)/10000])

subplot(312), title('Theta generator'), hold on, ylabel('Zp (Hz)'), xlabel('time (s)')
plot(t,zpt,'k','LineWidth',0.75)
xlim([0 length(t)/10000])

subplot(313), title('Layer 2'), hold on, ylabel('Zp (Hz)'), xlabel('time (s)')
plot(t,(zp2(ep1(1:end),:)),'b','LineWidth',1)
plot(t,(zp2(ep2(1:end),:)),'Color',[255 128 0]/255,'LineWidth',1)
plot(t,(zp2(ep3(1:end),:)),'g','LineWidth',1)
plot(t,(zp2(ep4(1:end),:)),'r','LineWidth',1)
plot(t,(zp2(ep5(1:end),:)),'k','LineWidth',1)

plot(t,(zp2(ep6(1:end),:)),'b--','LineWidth',0.75,'LineWidth',1)
plot(t,(zp2(ep7(1:end),:)),'--','Color',[255 128 0]/255,'LineWidth',1)
plot(t,(zp2(ep8(1:end),:)),'g--','LineWidth',0.75,'LineWidth',1)
plot(t,(zp2(ep9(1:end),:)),'r--','LineWidth',0.75,'LineWidth',1)
plot(t,(zp2(ep10(1:end),:)),'k--','LineWidth',0.75,'LineWidth',1)

plot(t,zp2(ep11(1:end),:),'b:','LineWidth',0.75,'LineWidth',1.5)
plot(t,zp2(ep12(1:end),:),':','Color',[255 128 0]/255,'LineWidth',1.5)
plot(t,zp2(ep13(1:end),:),'g:','LineWidth',0.75,'LineWidth',1.5)
plot(t,zp2(ep14(1:end),:),'r:','LineWidth',0.75,'LineWidth',1.5)
plot(t,zp2(ep15(1:end),:),'k:','LineWidth',0.75,'LineWidth',1.5)

legend('ep1','','','','ep2','','','','','','ep3','','','','','ep4','','','','','','ep5','','','','ep6','','','','ep7','','','','','','ep8','','','','','ep9','','','','','','ep10','','','','ep11','','','','ep12','','','','','','ep13','','','','','ep14','','','','','','ep15','','','')
xlim([0 length(t)/10000])


%
figure(2)

subplot(3,1,1), title('Layer 2'), hold on, ylabel('Zp (Hz)'),
plot(t,(zp2(ep2(1:2),:)),'m','LineWidth',1.5)
plot(t,(zp2(ep3(3),:)),'Color',[200 100 16]/256,'LineWidth',1.5)
plot(t,(zp2(ep9(4),:)),'c','LineWidth',1.5)
plot(t,(zp2(ep1,:)),'b','LineWidth',0.5)
plot(t,(zp2(ep2(3:end),:)),'Color',[255 128 0]/255,'LineWidth',0.5)
plot(t,(zp2(ep3([1:2 4:end]),:)),'g','LineWidth',0.5)
plot(t,(zp2(ep4,:)),'r','LineWidth',0.5)
plot(t,(zp2(ep5,:)),'k','LineWidth',0.5)
%new from 2nd sequence
plot(t,(zp2(ep6,:)),'b--','LineWidth',0.5)
plot(t,(zp2(ep7([1 3:end]),:)),'--','Color',[255 128 0]/255,'LineWidth',0.5)
plot(t,(zp2(ep8,:)),'g--','LineWidth',0.5)
plot(t,(zp2(ep9([1:3 5:end]),:)),'r--','LineWidth',0.5)
plot(t,(zp2(ep10,:)),'k--','LineWidth',0.5)
%new from 3rd sequence
plot(t,(zp2(ep11,:)),'b:','LineWidth',1)
plot(t,(zp2(ep12([1 3:end]),:)),':','Color',[255 128 0]/255,'LineWidth',1)
plot(t,(zp2(ep13([1:2 4:end]),:)),'g:','LineWidth',1)
plot(t,(zp2(ep14(1:4),:)),'r:','LineWidth',1)
plot(t,(zp2(ep15,:)),'k:','LineWidth',1)
legend('common Ep 2-14','','common Ep 3-7-12')
xlim([0.28 0.46])
ylim([0 5]), yticks([1:5])

subplot(3,1,2), hold on, ylabel('Zp (Hz)'),
plot(t,(zp2(ep3(3),:)),'Color',[200 100 16]/256,'LineWidth',1.5)
plot(t,(zp2(ep9(4),:)),'c','LineWidth',1.5)
plot(t,(zp2(ep2(1:2),:)),'m','LineWidth',1.5)
plot(t,(zp2(ep6,:)),'b--','LineWidth',0.5)
plot(t,(zp2(ep7([1 3:end]),:)),'--','Color',[255 128 0]/255,'LineWidth',0.5)
plot(t,(zp2(ep8,:)),'g--','LineWidth',0.5)
plot(t,(zp2(ep9([1:3 5:end]),:)),'r--','LineWidth',0.5)
plot(t,(zp2(ep10,:)),'k--','LineWidth',0.5)
%new from 1st sequence
plot(t,(zp2(ep1,:)),'b','LineWidth',0.5)
plot(t,(zp2(ep2(3:end),:)),'Color',[255 128 0]/255,'LineWidth',0.5)
plot(t,(zp2(ep3([1:2 4:end]),:)),'g','LineWidth',0.5)
plot(t,(zp2(ep4,:)),'r','LineWidth',0.5)
plot(t,(zp2(ep5,:)),'k','LineWidth',0.5)
%new from 3rd sequence
plot(t,(zp2(ep11,:)),'b:','LineWidth',1)
plot(t,(zp2(ep12([1 3:end]),:)),':','Color',[255 128 0]/255,'LineWidth',1)
plot(t,(zp2(ep13([1:2 4:end]),:)),'g:','LineWidth',1)
plot(t,(zp2(ep14(1:4),:)),'r:','LineWidth',1)
plot(t,(zp2(ep15,:)),'k:','LineWidth',1)
legend('common Ep 3-7-12','common Ep 9-13')
xlim([0.52 0.72])
ylim([0 5]), yticks([1:5])

subplot(3,1,3), title(''), hold on, ylabel('Zp (Hz)'), xlabel('time (s)')
plot(t,(zp2(ep3(3),:)),'Color',[200 100 16]/256,'LineWidth',1.5)
plot(t,(zp2(ep9(4),:)),'c','LineWidth',1.5)
plot(t,(zp2(ep2(1:2),:)),'m','LineWidth',1.5)
plot(t,(zp2(ep11,:)),'b:','LineWidth',1)
plot(t,(zp2(ep12([1 3:end]),:)),':','Color',[255 128 0]/255,'LineWidth',1)
plot(t,(zp2(ep13([1:2 4:end]),:)),'g:','LineWidth',1)
plot(t,(zp2(ep14(1:4),:)),'r:','LineWidth',1)
plot(t,(zp2(ep15,:)),'k:','LineWidth',1)
%new fron 2nd sequence
plot(t,(zp2(ep6,:)),'b--','LineWidth',0.5)
plot(t,(zp2(ep7([1 3:end]),:)),'--','Color',[255 128 0]/255,'LineWidth',0.5)
plot(t,(zp2(ep8,:)),'g--','LineWidth',0.5)
plot(t,(zp2(ep9([1:3 5:end]),:)),'r--','LineWidth',0.5)
plot(t,(zp2(ep10,:)),'k--','LineWidth',0.5)
%new from 1st sequence
plot(t,(zp2(ep1,:)),'b','LineWidth',0.5)
plot(t,(zp2(ep2(3:end),:)),'Color',[255 128 0]/255,'LineWidth',0.5)
plot(t,(zp2(ep3([1:2 4:end]),:)),'g','LineWidth',0.5)
plot(t,(zp2(ep4,:)),'r','LineWidth',0.5)
plot(t,(zp2(ep5,:)),'k','LineWidth',0.5)
legend('common Ep 3-7-12','common Ep 9-13','common Ep 2-14')
xlim([1.02 1.18])
ylim([0 5]), yticks([1:5])

