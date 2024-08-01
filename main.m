 close all;
clear all;
clc
clear
% Infrared small target detction based on independent component analysis on
% 5-D spatial-temporal domain
%% parameter setting
format long
addpath('utils/');
addpath('metric_utils\');
addpath('tensor_toolbox\');

readPath = '.\image';
savePath = '.\result';

if ~exist(savePath)
    mkdir(savePath);
end

frame = 7;
patchSize=20; 
lambdaL=15; 
mu = 0.0007;
x=0.007;

tuneopts.temporal_step = frame;
tuneopts.patchSize = patchSize;
tuneopts.lambdaL = lambdaL;
tuneopts.mu=mu; 
tuneopts.x=x; 
 
target_detection(char(readPath), savePath, tuneopts); 
