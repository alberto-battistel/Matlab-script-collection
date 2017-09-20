%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting of EIS data
% first with always convergent algorithm (Simulated Annealing SA), but slow
% and then refining with a simplex (SX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

%% Z are in Data
% collected into cells

cycleN = 1;

col_Z = [2,3];
col_F = 1;

% restrict frequencies
l = 11:35;

f = Data{cycleN}(l,col_F);
Z_exp = Data{cycleN}(l,col_Z(1)) - 1j*Data{cycleN}(l,col_Z(2));


% R1 = 2.1e-3;
% R2 = 5e-3;
% Q = 1e-1;
% a = 0.9;
% Z_exp = R1+(1./R2 + (2*1j*pi*f).^a*Q).^-1;


%% Equivalent circuit
% it is an external function
% the parameters to seek are in log10 scale to improve relative weight
fun = str2func('RsCPERctCPE');

%% Starting guess
% starting guess for Simulated Annealing SA
% they have to be the same number requested by the equivalent circuit
% function
V0_SA = [-3; -2; -1; -1; -1; -1];

%% Loss function
% the same for SA and SX
n = length(f)-length(V0_SA);
lossFun = @(V) sum( ((real(fun(V,f)) - real(Z_exp)).^2 + ...
        (imag(fun(V,f)) - imag(Z_exp)).^2) ./ (n*(real(Z_exp).^2 + (imag(Z_exp).^2))));

%% Optimization options
% options for simulated annealing
options_SA = anneal();

     options_SA.CoolSched = @(T)(.95*T);
     options_SA.Generator = @(x)(x(:)'+(randperm(length(x))==length(x))*randn/100)';
      options_SA.InitTemp = 1000;
    options_SA.MaxConsRej = 400;
    options_SA.MaxSuccess = 80;
      options_SA.MaxTries = 5000;
      options_SA.StopTemp = 1.0000e-15;
       options_SA.StopVal = 3e-3;
     options_SA.Verbosity = 1;

% options for simplex
options_SX = optimset('Display','none','TolFun',1e-12,'TolX',1e-12,'MaxFunEvals',1e5,'MaxIter',1e4);

%% minimization
% SA
[Vf_SA,chi2_SA] = anneal(lossFun, V0_SA, options_SA);

% SX
% start from the parameters found by SA
[Vf_SX,chi2_SX] = fminsearch(@(V) lossFun(V), Vf_SA, options_SX);

%% Z calculated
% SA
Z_SA = fun(Vf_SA,f);

% SX and parameters in linear scale
[Z_SX, V_lin] = fun(Vf_SX,f);

%% Plot
% Nyquist
figure(1)
clf
hold on

% experimental data
plot(conj(Z_exp), 'b.')

% SA
plot(conj(Z_SA), 'g-')

% SX
plot(conj(Z_SX), 'r-')

hold off
axis equal
xlabel('Z_{Re}')
ylabel('-Z_{Im}')
legend('Z_{Exp}', 'Z_{SA}', 'Z_{SX}')

% Bode
figure(2)
clf
subplot(2,1,1)
hold on 
loglog(f, abs(Z_exp), 'b.')
loglog(f, abs(Z_SA), 'g-')
loglog(f, abs(Z_SX), 'r-')
hold off
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('f')
ylabel('Abs(Z)')
legend('Z_{Exp}', 'Z_{SA}', 'Z_{SX}')

subplot(2,1,2)
hold on
Phase = @(Z) atan(imag(Z)./real(Z))/pi*180;
semilogx(f, -Phase(Z_exp), 'b.')
semilogx(f, -Phase(Z_SA), 'g-')
semilogx(f, -Phase(Z_SX), 'r-')
hold off
set(gca, 'XScale', 'log')
xlabel('f')
ylabel('-Phase(Z)')




