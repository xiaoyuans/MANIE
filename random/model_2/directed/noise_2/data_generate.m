clc;
clear;
close all;
addpath('Models/','Functions/')

MODEL={'kuramoto1','kuramoto2','michaelis_menten','roessler'};
BASIS={'polynomial','polynomial_diff','fourier','fourier_diff','power_series','RBF'};

N=25;
NI=4;
S=50;
M=5;

simulate2(MODEL{2},N,NI,S,M);