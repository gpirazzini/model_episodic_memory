%% Main

clear all
close all
clc

Npop = 75; % n° of populations

% delays within and between layer
D_intraLayer = 1;
D_extraLayer = 30;

TRAINING_L1_ortho          % comment for training ORTHOGONAL sequences
% TRAINING_L1_shared       % uncomment for training SHARED sequences

% figure summarizing the 4 types of trained synapses
figure()
subplot(221)
imagesc(Wp_L1L1), colorbar, title('Wp_L_1_L_1')
subplot(222)
imagesc(Wf_L1L1), colorbar, title('Wf_L_1_L_1')
subplot(223)
imagesc(Af_L1L1), colorbar, title('Af_L_1_L_1')
subplot(224)
imagesc(Wp_L1L2), colorbar, title('Wp_L_1_L_2')
sgtitle('Trained synapses')

save synapses_ORTHO    Wp_L1L1 Wf_L1L1 Af_L1L1 Wp_L1L2 %uncomment for save the results (ORTHO)
% save synapses_SHARED Wp_L1L1 Wf_L1L1 Af_L1L1 Wp_L1L2 %uncomment for save the results (SHARED)

