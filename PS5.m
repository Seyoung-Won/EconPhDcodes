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

%% Q1.2 - deterministic RA version


r = 1/beta - 1 ;







%% simulate AR1 process for the log labour supply log l

l0 = 0;
T  = 250;

xt1 = simulate_AR1(T,mu,phi,sigeps,l0);

[chain,state] = simulate_markov_chain(T,nodes,P,l0);

xt2 = chain;

figure(1)
plot([xt1,xt2])
legend('AR1','markov chain')
xlabel('time')
ylabel('x(t)')




return











