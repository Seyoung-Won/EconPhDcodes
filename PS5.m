clear ;
close all;
clc;

%%%%% AR1 process
% Parameters

beta = 0.96 ;
gamma = 1.5 ;
rho = 0.9 ;
alpha = 0.36 ;
delta = 0.10 ;

% epsilon is normally distributed with:
mu     = 0;
sigma = 0.2 ;

%%%% discrete-state approximation to AR1 

N = 7 ;         % Nodes of the grid = 7
m = 2.5 ;       % standard deviation to discretize AR(1) = 2.5

[logell, ell_prob] = tauchen(N,mu,rho,sigma,m) ;

ell = exp(logell) ;

[eigvec, eigval] = eig(ell_prob') ;
[~, arg] = min(abs(diag(eigval) - 1)) ;
unit_eigvec = eigvec(:, arg) ;

% Stationary distribution vector of the transition matrix
Zprob_bar = unit_eigvec / sum(unit_eigvec) ;

expected_ell = ell' * Zprob_bar ;

% Alternatively you can use: 
% P = mc_invdist(ell_prob) ; % which gives the same as Zprob_bar

%%%%%%%%%% Q1.2 - deterministic RA version

% at steady state (1+r)*beta = 1, we can find such r
% given this r value, we know that at steady state r + delta = MPK 
% which is the FOC condition of the capital dynamics: 
% interest rate to wealth r = MPK - depreciation 
% MPK = alpha* (K/L)^(alpha - 1) = r + delta, rearranging gives
% K = (alpha / (r + delta))^(1/(1 - alpha))*L , at steady state L = E(ell)

r_ss = 1/beta - 1 ;
K_RA = (alpha / (r_ss + delta))^(1/(1-alpha)) * expected_ell ;

%%%%%%%%%% Q1.3 - Calculating the General Equilibrium numerically

amin = 0 ;
amax = 50 ;
nbgrid = 150 ;
nbstate = 7 ;
agrid = linspace(amin, amax, nbgrid) ;
