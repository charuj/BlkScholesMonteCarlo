% Pricing a European option using Black-Scholes formula and Monte Carlo simulations

S0 = 100;     % price of underlying stock today
K = 105;      % strike at expiry
mu = 0.05;    % expected return
sigma = 0.2;  % volatility
r = 0.05;     % risk-free rate
T = 1.0;      % years to expiry
nruns = 10000; % number of simulated paths
dt = 1/365;   % time steps
etime = 365;   % days to expiry

% Implementing Black-Scholes pricing formula

[c, p]= blsprice(S0, K, r, T,sigma);  % I use the native Black-Scholes function in Matlab 

% the blsprice computes the premium (c,p) that then needs to be added to the original price (S0)

putBS_Price= p + S0
callBS_Price= c + K*exp(-r*T)

% Implementing Monte Carlo pricing procedure

%[callMC_Price, putMC_Price] = MC_price(S0, K, T, r, mu, sigma);


final_prices_list=zeros(nruns,1);

for i=1:nruns
    Wpast= 0;
    S_current=S0;
    S_current_history=zeros(etime,1);
    for calculations= 1:etime;
       W_current= normrnd(0,1);
       delta_W= W_current-Wpast;
     
       sig_sq= sigma^2 ;
       sig_sq=sig_sq/2;
       arg1= (mu-sig_sq)*dt;
       root_dt= sqrt(dt);
       arg2= sigma*root_dt*delta_W;
     
       S_current=S_current*exp(arg1 + arg2);
       S_current_history(calculations)=S_current;
       
       Wpast= W_current ;
    end
    
    plot(S_current_history);  % PLOTS FIGURE 1 
    hold on;
    final_prices_list(i)=S_current;
end 
%S = GRWPaths(S0, mu, sig, T, etime, nruns);

% Calculate the payoff for each path for a Put
%PutPayoffT = max(K-mean(S),0);
PutPayoffT = max(K-(final_prices_list),0)


% Calculate the payoff for each path for a Call
CallPayoffT = max((final_prices_list)-K,0)

% Discount back
putPrice = mean(PutPayoffT)*exp(-r*T);
callPrice = mean(CallPayoffT)*exp(-r*T);

% add the discounts to the original price (S0)
putPrice= putPrice + S0
callPrice= callPrice + K*exp(-r*T)

