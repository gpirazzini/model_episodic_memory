%% MAIN

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

% time parameters
t_sim=9;
dt=0.0001;
t=0:dt:t_sim;
T=length(t);

% load trained synapses
% load <<INSERT YOUR DIRECTORY HERE>> ... synapses_SHARED


% IMAGINATION: *3/3 ALL synapses ...depending on the operation mode
% DREAMING: *2/3 ALL synapses
% SCHIZOPHRENIA: *1.5/3 only Af_L1L1 synapses
Wp_L1L1 = Wp_L1L1*3/3;
Wf_L1L1 = Wf_L1L1*3/3;
Af_L1L1 = Af_L1L1*3/3;
Wp_L1L2 = Wp_L1L2*3/3;
Wp_L2L1=eye(Npop)*300*3/3;
 
% IMAGINATION: 150+150*rand
% DREAMING: 250+250*rand
% SCHIZOPHRENIA: 250+250*rand
mediaIN=150+150*rand([Npop,length(t)]); %noise 

ISOLATION_sim


%%  Figures

figure(1)

subplot(311), title('Layer 2'), hold on, ylabel('Zp (Hz)'), ylim([0 5])
plot(t,(sum(zp2(ep1(1:end),:))/(size(ep1,2))),'b')
plot(t,(sum(zp2(ep2(1:end),:))/(size(ep2,2))),'Color',[255 128 0]/255),
plot(t,(sum(zp2(ep3(1:end),:))/(size(ep3,2))),'g')
plot(t,(sum(zp2(ep4(1:end),:))/(size(ep4,2))),'r')
plot(t,(sum(zp2(ep5(1:end),:))/(size(ep5,2))),'k')
plot(t,(sum(zp2(ep6(1:end),:))/(size(ep6,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep7(1:end),:))/(size(ep7,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep8(1:end),:))/(size(ep8,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep9(1:end),:))/(size(ep9,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep10(1:end),:))/(size(ep10,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep11(1:end),:))/(size(ep11,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep12(1:end),:))/(size(ep12,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep13(1:end),:))/(size(ep13,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep14(1:end),:))/(size(ep14,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep15(1:end),:))/(size(ep15,2))),':','LineWidth',1.5)
xlim([0 3])

subplot(312), hold on, ylabel('Zp (Hz)'), ylim([0 5])
plot(t,(sum(zp2(ep1(1:end),:))/(size(ep1,2))),'b')
plot(t,(sum(zp2(ep2(1:end),:))/(size(ep2,2))),'Color',[255 128 0]/255),
plot(t,(sum(zp2(ep3(1:end),:))/(size(ep3,2))),'g')
plot(t,(sum(zp2(ep4(1:end),:))/(size(ep4,2))),'r')
plot(t,(sum(zp2(ep5(1:end),:))/(size(ep5,2))),'k')
plot(t,(sum(zp2(ep6(1:end),:))/(size(ep6,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep7(1:end),:))/(size(ep7,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep8(1:end),:))/(size(ep8,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep9(1:end),:))/(size(ep9,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep10(1:end),:))/(size(ep10,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep11(1:end),:))/(size(ep11,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep12(1:end),:))/(size(ep12,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep13(1:end),:))/(size(ep13,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep14(1:end),:))/(size(ep14,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep15(1:end),:))/(size(ep15,2))),':','LineWidth',1.5)
xlim([3 6])

subplot(313), hold on, ylabel('Zp (Hz)'), xlabel('time (s)'), ylim([0 5])
plot(t,(sum(zp2(ep1(1:end),:))/(size(ep1,2))),'b')
plot(t,(sum(zp2(ep2(1:end),:))/(size(ep2,2))),'Color',[255 128 0]/255),
plot(t,(sum(zp2(ep3(1:end),:))/(size(ep3,2))),'g')
plot(t,(sum(zp2(ep4(1:end),:))/(size(ep4,2))),'r')
plot(t,(sum(zp2(ep5(1:end),:))/(size(ep5,2))),'k')
plot(t,(sum(zp2(ep6(1:end),:))/(size(ep6,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep7(1:end),:))/(size(ep7,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep8(1:end),:))/(size(ep8,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep9(1:end),:))/(size(ep9,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep10(1:end),:))/(size(ep10,2))),'--','LineWidth',1.25)
plot(t,(sum(zp2(ep11(1:end),:))/(size(ep11,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep12(1:end),:))/(size(ep12,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep13(1:end),:))/(size(ep13,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep14(1:end),:))/(size(ep14,2))),':','LineWidth',1.5)
plot(t,(sum(zp2(ep15(1:end),:))/(size(ep15,2))),':','LineWidth',1.5)
legend('ep1','ep2','ep3','ep4','ep5','ep6','ep7','ep8','ep9','ep10','ep11','ep12','ep13','ep14','ep15')
xlim([6 9])

