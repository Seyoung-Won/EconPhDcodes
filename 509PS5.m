% Macro 509 Problem Set 5 - Seyoung Won
% This code definitely WORKS but it may take you 10 minutes to get the result.
clearvars ;
clc ;

%%%%% AR1 process Parameters
beta = 0.96 ;
gamma = 1.5 ;
rho = 0.9 ;
alpha = 0.36 ;
delta = 0.10 ;

% error is normally distributed with:
mu     = 0 ;
sigma = 0.2 ;

N = 7 ;         % Nodes of the grid = 7
m = 2.5 ;       % standard deviation to discretize AR(1) = 2.5

% Use tauchen method to discritize the labour into grids
% It gives N=7 gridpoints of log labour and labour transition matrix 
[logell, ell_prob] = tauchen(N,mu,rho,sigma,m) ;    
ellgrid = exp(logell) ;

ell_stdist = mc_invdist(ell_prob) ;     % Stationary Distribution of labour
expected_ell = ellgrid' * ell_stdist ;  % Expectation of labour

%%%%%%%%%%%%%%%%%%%% Q1.2 - deterministic RA version %%%%%%%%%%%%%%%%%%%%

% at steady state (1+r)*beta = 1, we can find such r
% given this r value, we know that at steady state r + delta = MPK 
% which is the FOC condition of the capital dynamics: 
% interest rate to wealth r = MPK - depreciation 
% MPK = alpha* (K/L)^(alpha - 1) = r, rearranging gives
% K = (alpha / MPK)^(1/(1 - alpha))*L , at steady state L = E(ell)

r = 1/beta - 1 + delta ;         % this is MPK
K_RA = (alpha / r)^(1/(1 - alpha))* expected_ell ;
w = (1-alpha)*expected_ell^(-alpha)*K_RA^(alpha) ;

%%%%%%%%%% Q1.3 - Calculating the General Equilibrium numerically

% Define Asset grids: amin=0, amax=50, no.grid=150
agrid = linspace(0, 50, 150) ;

% wealth a, with 7 states and 150 gridpoints
% labour ell, with 7 states and 150 gridpoints
[ell, a] = meshgrid(ellgrid, agrid) ;  

% dim. of all states = 7 states of labour from tauchen * 150 wealth grids
s_dim = max(size(ellgrid)) * max(size(agrid)) ;

% Markov chain: Vectors of state variables; 
% Note that in Prof. Greaney's notation s is NM by 1 but 1 is a tuple
% Hence we actually need a NM x 2 (a, l) vectors
s = [reshape(a, [s_dim, 1]), reshape(ell, [s_dim, 1])];

epsilon = 10^(-4) ;

K = K_RA;
crit = 1;
while crit > 0.001
    
% make null matrices : to be filled in by iterations
Vfn_old = zeros(max(size(agrid)), max(size(ellgrid))) ;
Vfn_new = zeros(max(size(agrid)), max(size(ellgrid))) ;
Policyfn = zeros(max(size(agrid)), max(size(ellgrid))) ;

max_diff = 1;
while max_diff > (1-beta)*epsilon
    for i = 1 : s_dim

        % calculatie consumption
        a_i = a(i) ;
        ell_j = ell(i) ;
        c = (1 + r )*a_i + w * ell_j - agrid ;

        % get the non-negative consumption only
        posi_c = find(c >= 0) ;

        % calculate utility
        util = c(posi_c).^(1-gamma)/(1-gamma) ;

        % get the expected value function
        ell_location = find(ell_j == ellgrid) ;
        EVfn = ell_prob(ell_location, :)*Vfn_old(posi_c, :)' ;
         
        % Get the value function and fill into Vfn_new
        Vfn_all = util + beta .* EVfn ;
        Vfn_max = max(Vfn_all) ;

        % update the value function
        Vfn_new(i) = Vfn_max ;

        % update the policy function
        Policyfn(i) = agrid(posi_c(Vfn_all == Vfn_max)) ;
    end
    
    max_diff = max(abs(Vfn_new - Vfn_old), [], 'all');
    fprintf('Maximum Difference : %i\n', max_diff) ;   
    Vfn_old = Vfn_new ; % updating the value function
    
end

    % Getting the state transition matrix using Markov Chain vector of states
    transition_matrix = zeros(s_dim, s_dim) ;
    for i = 1:s_dim
       for j = 1:s_dim 
           % get the transition probability, note s(:,2) = ell
           i_dim = find(s(i,2) == ellgrid) ; 
           j_dim = find(s(j,2) == ellgrid) ; 
           s_prob = ell_prob(i_dim, j_dim) ; % state transition probability

           % Get the policy function g(s) = g(i(s), j(s))
           % i(s) = mod(s-1, M) +1 ; j(s) = floor((s-1)/M) + 1 ; M = 150
           % Probability "exists" if g(s) = i(s')
           g_s = Policyfn(mod(i-1,150) + 1, floor((i-1)/150) + 1) ; % chosen level
           i_s = s(j,1) ;                                             % actual level

           % Put into the state transition matrix
           transition_matrix(i,j) = s_prob * (g_s == i_s) ;
       end
    end     

    % Get Sn the Aggregate Savings Sn = sum sum a_i * psi(a_i, l_i)
    % First get psi(a,l) = stationary distribution of the transition mat
    stationary_distribution = mc_invdist(transition_matrix) ;
    % Aggregate saving is a dot product of policyfn_asset * st.dist above
    Savings = sum(Policyfn .* reshape(stationary_distribution, 150, 7), 'all') ;

    % compute convergence criterion and update K, w=MPL, r=MPK
    crit = abs((K - Savings) / K) ; 
    K = K + 0.1*(Savings - K) ; 
    w = (1 - alpha) *expected_ell^(-alpha) * K^(alpha) ;
    r = alpha*expected_ell^(1-alpha)*K^(alpha-1) - delta ;

end

%
clc
fprintf('Q1.1: Expected labour = %f\n', expected_ell);
fprintf('Q1.2: Equilibrium Capital stock in RA model = %f\n', K_RA);
fprintf('Q1.3a: Aiyagari Capital stock = %f\n', K);
fprintf('Q1.3a: Aiyagari Wage = %f\n', w);
fprintf('Q1.3a: Aiyagari Interest Rate = %f\n', r);


% Plot the Value functions iterations
% Each of Value, Policy and stationary distributions correspond to 
% 7 eigenvalues of the labour grid from tauchen
figure;
hold on
legend()
for val = 1:7
    plot(agrid, Vfn_new(:, val),'DisplayName',string(ellgrid(val)));
end
xlabel('a')
ylabel('V(a,l)')
title('Value Functions Iteration')
legend('Location', 'Southeast')
hold off;
saveas(gcf,'value_functions.png');

% Plot the policy functions iterations
figure;
hold on
legend()
for val = 1:7
    plot(agrid, Policyfn(:, val),'DisplayName',string(ellgrid(val)));
end
xlim([0 50])
ylim([0 51])
xlabel('a')
ylabel('a''(a,l)') 
title('Policy Functions Iterations')
legend('Location', 'Southeast')
hold off;
saveas(gcf,'policy_functions.png');

% Plot density of state variables
density = reshape(stationary_distribution, 150, 7);
figure;
hold on
legend()
for val = 1:7
    plot(agrid, density(:, val),'DisplayName',string(ellgrid(val)));
end
xlabel('a')
ylabel('psi(a,l)')
title('Density of state functions - Stationary Distributions')
legend('Location', 'Northeast')
hold off;
saveas(gcf,'density_stationary_dist.png');
