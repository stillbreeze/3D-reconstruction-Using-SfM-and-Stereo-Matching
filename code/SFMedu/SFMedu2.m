clc;
disp('SFMedu: Structrue From Motion for Education Purpose');
disp('Version 2 @ 2014');
disp('Written by Jianxiong Xiao (MIT License).. Modifications by nick rhinehart');


%% set up things
clear;
close all;
addpath(genpath('matchSIFT'));
addpath(genpath('denseMatch'));
addpath(genpath('RtToolbox'));

%% parameters for visualization, allowing for focal length adjustment
% and graph merging strategy. change data_seq_idx to one of {0, 1, 2}
visualize = true;
adjust_focal_length = true;
do_sequential_match = false;
do_dense = true;
data_seq_idx = [0 3 4];
maxSize = 640;

% fprintf('running on head\n');
full_pipeline(data_seq_idx(1), adjust_focal_length, do_sequential_match, maxSize, visualize, do_dense);
% fprintf('running on iron\n');
full_pipeline(data_seq_idx(2), adjust_focal_length, do_sequential_match, maxSize, visualize, do_dense);
% fprintf('running on RI robo\n');
full_pipeline(data_seq_idx(3), adjust_focal_length, do_sequential_match, maxSize, visualize, do_dense);

% fprintf('running on head\n');
% full_pipeline(data_seq_idx, adjust_focal_length, do_sequential_match, maxSize, visualize, do_dense);
% 
% data_seq_idx = data_seq_idx + 1;
% fprintf('running on couch\n');
% full_pipeline(data_seq_idx, adjust_focal_length, do_sequential_match, maxSize, visualize, do_dense);
% 
% data_seq_idx = data_seq_idx + 1;
% fprintf('running on mug\n');
% full_pipeline(data_seq_idx, adjust_focal_length, do_sequential_match, maxSize, visualize, do_dense);