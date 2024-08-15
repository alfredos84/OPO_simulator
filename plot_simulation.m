clc; clear all;
% close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultAxesFontSize',20)


% sr = load ('output_signal_9999_r.dat'); si = load ('output_signal_9999_i.dat');
% ir = load ('output_idler_9999_r.dat'); ii = load ('output_idler_9999_i.dat');
% pr = load ('output_pump_9999_r.dat'); pi = load ('output_pump_9999_i.dat');

sr = load ('output_signal_r.dat'); si = load ('output_signal_i.dat');
ir = load ('output_idler_r.dat'); ii = load ('output_idler_i.dat');
pr = load ('output_pump_r.dat'); pi = load ('output_pump_i.dat');

sig = sr + 1i* si;  S = abs(sig).^2;
idl = ir + 1i* ii;  I = abs(idl).^2;
pump = pr + 1i* pi; P = abs(pump).^2;

figure(); hold on;
plot(P); plot(S); plot(I)