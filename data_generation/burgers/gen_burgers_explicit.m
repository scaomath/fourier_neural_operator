%% params
% u0 = GRF1(s/2, 0, gamma, tau, sigma, "periodic");
% s = 1024
clear;
m = 0;
n_space = 1024;
n_freq = n_space/2;
steps = 200;
gamma = 2.5;
tau = 7;
sigma = 7^(2);
visc = 1e-1;

const = 2*pi;
my_eigs = sqrt(2)*(abs(sigma).*((const.*(1:n_freq)').^2 + tau^2).^(-gamma/2));

xi_alpha = randn(n_freq,1);
alpha = my_eigs.*xi_alpha;

xi_beta = randn(n_freq,1);
beta = my_eigs.*xi_beta;
a = alpha/2;
b = -beta/2;

c = [flipud(a) - flipud(b).*1i; m + 0*1i;a + b.*1i];
%% make chebfun object
uu0 = chebfun(c, [0 1], 'trig', 'coeffs');
u0 = chebfun(@(t) uu0(t - 0.5), [0 1], 'trig');

%% u = burgers1(u0, tspan, s, visc)
tspan = linspace(0,1,steps+1);
S = spinop([0 1], tspan);
S.lin = @(u) + visc*diff(u,2);
S.nonlin = @(u) - 0.5*diff(u.^2);
S.init = u0;
u = spin(S,n_space,1e-4,'plot','on'); 

%%
output = zeros(1, steps, n_space);
for k=2:(steps+1)
    output(1,k,:) = u{k}.values;
end
plot(squeeze(output(1,end,:)))
set(gcf,'position',[100,500,800,500])
%%
% for i = 1:20:201
%     plot(squeeze(output(1,i,:)))
%     hold on
% end
waterfall(u)