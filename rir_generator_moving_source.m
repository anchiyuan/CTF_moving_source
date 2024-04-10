clc; clear;
close all;

%% RIR parameter %%
SorNum = 1;                                              % source number
MicNum = 30;                                             % number of microphone
c = 343;                                                 % Sound velocity (m/s)
fs = 16000;                                              % Sample frequency (samples/s)

% ULA %
MicStart = [1, 1.5, 1];
spacing = 0.02;
MicPos = zeros(MicNum, 3);
for i = 1:MicNum
    MicPos(i, :) = [MicStart(1, 1) + (i-1)*spacing, MicStart(1, 2), MicStart(1, 3)];
end

room_dim = [5, 6, 2.5];                                  % Room dimensions [x y z] (m)
SorPos = [2, 2.6, 1];                                    % first source position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SorPosNum = 4;                                           % number of source positions
Sor_spacing = 1;                                       % source spacing
reverberation_time = 0.6;                                % Reverberation time (s)
points_rir = 12288;                                      % Number of rir points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtype = 'omnidirectional';                               % Type of microphone
order = -1;                                              % -1 equals maximum reflection order!
dim = 3;                                                 % Room dimension
orientation = 0;                                         % Microphone orientation (rad)
hp_filter = 1;                                           % Disable high-pass filter

% calculate source position %
for j = 2:SorPosNum
    SorPos(j, :) = [SorPos(j-1, 1), SorPos(j-1, 2) + Sor_spacing, SorPos(j-1, 3)];
end

% 畫空間圖 %
figure(1);
plot3( [0 room_dim(1, 1) room_dim(1, 1) 0 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1) 0 0 room_dim(1, 1) room_dim(1, 1)], ...
       [0 0 room_dim(1, 2) room_dim(1, 2) 0 0 0 room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) room_dim(1, 2) 0 0 0], ...
       [0 0 0 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0 0 room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) room_dim(1, 3) 0] , 'k')
hold on
plot3(MicPos(:, 1), MicPos(:, 2), MicPos(:, 3), 'r.', 'MarkerSize', 10)
hold on
plot3(SorPos(:, 1), SorPos(:, 2), SorPos(:, 3), '*', 'MarkerSize', 20)
hold off
xlabel('x\_axis')
ylabel('y\_axis')
zlabel('z\_axis')
title('空間圖')
shg

%% generate ground-truth RIR (h) %%
% 產生 RIR 和存.mat 檔 %
for j = 1:SorPosNum
    h = rir_generator(c, fs, MicPos, SorPos(j, :), room_dim, reverberation_time, points_rir, mtype, order, dim, orientation, hp_filter);
    rir_filename_str = ['h\h_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), '.', string(j), '.mat'];
    rir_filemane = join(rir_filename_str, '');
    save(rir_filemane, 'h')
end

% 存 SorPos .mat 檔 %
SorPos_filename_str = ['h\SorPos_', string(reverberation_time), 'x', string(points_rir), 'x', string(Sor_spacing), 'x', string(SorPosNum), '.mat'];
SorPos_filemane = join(SorPos_filename_str, '');
save(SorPos_filemane, 'SorPos')
