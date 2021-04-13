// commitment vs. discretion in the new Keynesian model, code by Seyoung Won

// x = output gap, pi = inflation, i = interest rate, w = welfare
// e_pi = NKPC (supply) shock, e_x NIS (demand) shock
var x, pi, i, w, e_pi, e_x ;

// real AR(1) shocks to e_pi and e_x
varexo eta_pi eta_x ;

parameters beta, gamma, lambda, sigma, rho_pi, sigma_pi, rho_x, sigma_x ;

beta    = 0.96 ; // discount factor
gamma   = 0.10 ; // output gap weight
lambda  = 0.25 ; // NKPC coefficient
sigma   = 2 ;    // CARA coefficient

// AR(1) shocks

rho_pi      = 0.95 ;    // AR(1) coefficient for NKPC shock
rho_x       = 0.95 ;    // AR(1) coefficient for NIS shock
sigma_pi    = 0.01 ;    // Std. Dev. of eta_pi
sigma_x     = 0.005 ;   // Std. Dev. of eta_x

// New Keynesian Optimal Monetary Model starts here
model;

pi = beta*pi(+1) + lambda*x + e_pi ;        // NKPC
x = x(+1) - 1/sigma*(i - pi(+1)) + e_x ;    // NIS


// so far so good, under requires revision


// Five monetary policy 'targeting rule' : 1, 2, 3a, 3b, 3c
// 1 - Discretion, 2 - Commitment, 
// 3a - Discrete inflation targeting, pi = 0, 
// 3b - Discrete output gap targeting, x = 0, 
// 3c = Taylor Rule: i = phi_pi*pi + phi_x*x

// Remove '//' in front for activation

//                                     // 1 - Discretion
//                                     // 2 - Commitment
//                                     // 3 - Discrete inflation targeting, pi = 0
//                                     // 4 - Discrete output gap targeting, x = 0 
// i = 1.5*pi + 0.5*x                  // 5 - Taylor Rule with phi_pi = 1.5, phi_x = 0.5


// inflation identity equation?


// above requires revision

// law of motion for shocks 
e_pi = rho_pi * e_pi(-1) + eta_pi ;
e_x = rho_x * e_x(-1) + eta_x ;

end;
// Hereby New Keynesian optimal Monetary model ends

// Guess Steady State
initval; // Initial value for guess

x  = 0 ;
pi = 0 ;
i = 0 ; 
w = 0 ;
e_pi  = 0 ;
e_x = 0 ;

end;

// Declare the variance of the exogeneous shocks
shocks;

var e_pi = sigma_pi^2 ;
var e_x = sigma_x^2 ;

end;

// Get steady state results
steady ;

// one std shock to inflation / output and evolution for 40 periods
stoch_simul(order=1,irf=40);
