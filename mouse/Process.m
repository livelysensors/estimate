%% Initialize

load('combined.mat');
valff = [0 3.4 5.4 7.8 10.6 15.0 20.1 24.4 30.4 40.1 50.6 100];

selff = [1,2,3,4,5,6,7,8,9,10,11,12];
selff = [1,2,3,4,6,7,8,9,10,11,12];
capff = ['PDFF =  0.0';
         'PDFF =  3.4';
         'PDFF =  5.4';
         'PDFF =  7.8';
%         'PDFF = 10.6';
         'PDFF = 15.0';
         'PDFF = 20.1';
         'PDFF = 24.4';
         'PDFF = 30.4';
         'PDFF = 40.1';
         'PDFF = 50.6';
         'PDFF = 100.'];
     
% selff = [1,4,6,8,11,12];
% capff = ['PDFF =  0.0';
%          'PDFF =  7.8';
%          'PDFF = 15.0';
%          'PDFF = 24.4';
%          'PDFF = 50.6';
%          'PDFF = 100.'];

% selff = [1,8,11,12];
% capff = ['PDFF =  0.0';
%          'PDFF = 24.4';
%          'PDFF = 50.6';
%          'PDFF = 100.'];

% Number of FF samples being considered
dat   = s(:,selff);
valff = valff(selff);
numff = size(dat,2);
ffrowmax = ceil(sqrt(numff));
ffcolmax = ceil(numff/ffrowmax);

% Drop 1st echo, Trim to just 1ms
trim = 2:4096;
%trim = 2:1000;

temax = te(max(trim(:)));
dat = dat(trim,:);
te = te(trim,:);

%% Plot Data

figure(1);
plot(te,real(dat),'LineWidth',3);
legend(capff);
xlabel('TE');
ylabel('Signal (a.u.)');
axis([0 100 0 2000]);

%% Plot SemiLog (monoexponential)

figure(2);
semilogy(te,real(dat),'LineWidth',3');
legend(capff);
xlabel('TE');
ylabel('Signal (a.u.)');
axis([0 100 10^1 10^3.5]);

%% Fitting DAEUN
% 
% figure(4);
% t2_bin = logspace(log10(10),log10(1000),1000);
% lambda = 1e-10;
% [f, y_fit] = ILT_1D_regul(abs(dat(:,4)),te,lambda,t2_bin);
% 
% 
%% Fitting MonoExponential Krishna

figure(5);
clear parms;
for sig=1:size(dat,2)
    
    moexp = fit(te,real(dat(:,sig)),'exp1','StartPoint',[2000 -1/100]);
    subplot(ffrowmax,ffcolmax,sig);
    plot(moexp,te,real(dat(:,sig)));
    axis([0 temax 0 2000]);
    xlabel('TE'); ylabel('Signal'); title(capff(sig,:));
    
    parms(sig,:) = [-1/moexp.b moexp.a];
end

figure(6);
% stem(parms(:,1),parms(:,2),'LineWidth',3);
% axis([0 250 0 2000]);
% xlabel('T2 (ms)');
% ylabel('Signal (a.u.)');
% title('Monoexponential Fit');
% legend(capff);
stem(valff,parms(:,1),'LineWidth',3);
xlabel('Fat Fraction (%)');
ylabel('Apparent T2 (ms)');
title('Monoexponential Fit');


%% Fitting BiExponential Krishna
% 
% figure(7);
% for sig=1:size(dat,2)
%     
%     biexp = fit(te,real(dat(:,sig)),'exp2','StartPoint',[1000 -1/100 500 -1/50]);
%     subplot(ffrowmax,ffcolmax,sig);
%     plot(biexp,te,real(dat(:,sig)));
%     axis([0 temax 0 2000]);
%     xlabel('TE'); ylabel('Signal'); title(capff(sig,:));
%     
%     parameters(sig,:) = [-1/biexp.b -1/biexp.d biexp.a biexp.c ];
% end
% 
% figure(8);
% stem(parameters(:,1:2).',parameters(:,3:4).','LineWidth',3);
% axis([0 250 0 2000]);
% xlabel('T2 (ms)');
% ylabel('Signal (a.u.)');
% title('Biexponential Fit');
% legend(capff);



