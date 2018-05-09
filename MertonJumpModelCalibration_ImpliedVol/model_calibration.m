clc
clear
close all

data = xlsread('data.xlsx'); % input data from Excel
S0 = data(1,1); % initial value of underlying
r = data(5,1); % interest rate
q = data(3,1); % dividend yield
maturity = data(:,2); % all option maturities
strike = data(:,3); % all option strikes
Pobs = data(:,4); % all market option prices

sigma_0 = 0.2; lambda_0 = 0.7; mu_j_0 = 0.0; sigma_j_0 = 0.2; % theta_0: initial guesses of parameter values
theta_0 = [sigma_0 lambda_0 mu_j_0 sigma_j_0]; % vector of parameter values. theta_0(1) = sigma_0, theta_0(2) = lambda_0, theta_0(3) = mu_j_0, theta_0(4) = sigma_j_0.
lb = [0  0  -5  0]; % parameter lower bounds
ub = [10  10  5  10]; % parameter upper bounds

opts = optimset('Display','iter','MaxFunEvals',2000,'TolFun',1E-5); % options for numerical solution minimization problem

T = unique(maturity);
NT = length(T);
theta_star = zeros(NT,length(theta_0));
figure
for m = 1:NT
    K = strike(maturity == T(m));
    
    theta_star(m,:) = lsqnonlin(@(theta) ( P(T(m), K, r, q, S0, theta(1), theta(2), theta(3), theta(4))' - Pobs(maturity == T(m)) ),...
        theta_0, lb, ub, opts); % opti
    mal model parameter set
    
    P_theta_star = P(T(m), K, r, q, S0, theta_star(m,1), theta_star(m,2), theta_star(m,3), theta_star(m,4)); % optimal model option prices
    
    P_theta_star_IV = blsimpv(S0, K, r, T(m), P_theta_star', [], q, [], {'put'}); % implied volatilities from optimal model option prices
    
    Pobs_IV = blsimpv(S0, K, r, T(m), Pobs(maturity == T(m)), [], q, [], {'put'}); % implied volatilities from market option prices

    subplot(1,NT,m)
    plot(K, P_theta_star_IV, '-s', K, Pobs_IV, '-o')
    title(['Model-implied vol. vs. Market-implied vol., $T = \left. {}\right.$', num2str(round(T(m),2))],'interpreter','latex')
    legend('Model-implied vol.','Market-implied vol.')
    ylabel('Implied volatility','interpreter','latex')
    xlabel('$K$','interpreter','latex')
end

disp('MJD optimal parameter values for different maturities')
disp('         T    sigma    lambda    mu_j    sigma_j')
disp([T theta_star])

function output = P(T, K, r, q, S0, sigma, lambda, mu_j, sigma_j)

N = 2^13; % number of terms in the cosine expansion
CallPut = -1; % = 1 for call; -1 for put
a = -10;
b = max(log(K))+10;

NK = length(K);
output = zeros(1,NK);

for l = 1:NK
    if CallPut == 1 % call option
        c = log(K(l)); d = b;
    elseif CallPut == -1 % put option
        c = a; d = log(K(l));
    end
    
    % cosine series coefficients
    chi = (1./(1+((0:N-1)*pi/(b-a)).^2)).*(cos((0:N-1)*pi*(d-a)/(b-a))*exp(d)-...
        cos((0:N-1)*pi*(c-a)/(b-a))*exp(c)+((0:N-1)*pi/(b-a)).*sin((0:N-1)*pi*(d-a)/(b-a))*exp(d)-...
        ((0:N-1)*pi/(b-a)).*sin((0:N-1)*pi*(c-a)/(b-a))*exp(c));
    psi(1) = d-c;
    psi(2:N) = (sin((1:N-1)*pi*(d-a)/(b-a))-sin((1:N-1)*pi*(c-a)/(b-a)))*(b-a)./((1:N-1)*pi);
    
    % payoff series coefficient
    if CallPut == 1 % call
        V = (2/(b-a))*(chi-K(l)*psi);
    elseif CallPut == -1 % put
        V = -(2/(b-a))*(chi-K(l)*psi);
    end
    
    % characteristic function of the Merton jump diffusion (MJD)
    phi = MJD_phi((0:N-1)*pi/(b-a), T, r, q, S0, sigma, lambda, mu_j, sigma_j);
    
    output(l) = exp(-r*T)*sum(real(exp(-1i*(0:N-1)*pi*a/(b-a)).*phi).*V.*([0.5 ones(1,N-1)]));
end
end

% MJD characteristic function
function output = MJD_phi(u, t, r, q, S0, sigma, lambda, mu_j, sigma_j)
g = r-q-sigma^2/2-lambda*(exp(mu_j+sigma_j^2/2)-1);
output = S0.^(1i*u).*exp(1i*g*t*u+(-u.^2*sigma^2/2+lambda*(exp(1i*mu_j*u-u.^2*sigma_j^2/2)-1))*t);
end