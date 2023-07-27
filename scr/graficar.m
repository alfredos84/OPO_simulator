clc; clear all; close all;
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');  
set(0,'defaultAxesFontSize',20)


folder = 'Simulations/cw_3eqs_PPLN_beta_0.8_N_4_GDD_100_LP_532nm';
T      = load([folder,'/T.dat']);
trt    = T(end)-T(1); % ps
Tp     = load([folder,'/Tp.dat']);
Tp     = Tp + max(Tp);
f      = load([folder,'/freq.dat']);

signal_r=load([folder,'/signal_output_r.dat']);
signal_i=load([folder,'/signal_output_i.dat']);

pump_r=load([folder,'/pump_output_r.dat']);
pump_i=load([folder,'/pump_output_i.dat']);


SIGNAL  = signal_r + 1j*signal_i;
PUMP   = pump_r  + 1j*pump_i;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C    = 299792458*1e6/1e12; % speed of ligth in vacuum [um/ps]
EPS0 = 8.8541878128e-12*1e12/1e6; %vacuum pertivity [W.ps/V²μm] 
np = 2.22515;    ns = 2.14883;     ni= ns;

waist       = 55; % beam waist radius [um]
spot        = pi*waist^2; % spot area [μm²]

cp = .5 * EPS0 * C * spot * np;
cs = .5 * EPS0 * C * spot * ns * sqrt(2);

Ps    = cs*abs(SIGNAL(end-length(T)+1:end)).^2;
Pp    = cp*abs(PUMP(end-length(T)+1:end)).^2;


idler = false;
if(idler)
    idler_r = load([folder,'/idler_output_r.dat']);
    idler_i = load([folder,'/idler_output_i.dat']); 
    IDLER  = idler_r + 1j*idler_i;
    ci = .5 * EPS0 * C * spot * ni * sqrt(2);
    Pi    = ci*abs(IDLER(end-length(T)+1:end)).^2;
end



h = figure('units','normalized','outerposition',[0 0 1 1]);    
hold on
plot( T, Ps )
if(idler)
    plot( T, Pi )
end
ylabel('Output signal power (W)')
xlabel('T (ps)')
ax= gca; ax.PlotBoxAspectRatio = [2,1,1];
box on; grid on;