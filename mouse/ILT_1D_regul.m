function [f, y_fit] = ILT_1D_regul(y, te, lambda, t2);


% This code is implemnetd to solve the following equation 
% using the non-negative least squares (NNLS) algorithm
% ||y - Kf ||_2^2 + lambda ||Rf||_2^2 s.t. f>=0
% R: Regularization matrix to impose spectral smoothness

%Input
% y: measured cpmg signal (n x 1)
% te: echo time (n x 1)
% lambda: regularization parameter for spectral smoothness
% (defaul: 1e-10)
% t2: dictionary parameters
% (defalut: t2 = logspace(log10(2.77), log10(300), 300); )
% 
%Ouput
%f: 1D spectrum of t2
%

if(nargin < 4)
    t2 = logspace(log10(2.77), log10(300), 300);
end
if(nargin < 3)
    lambda = 1e-10;
end

y = y(:);
te = te(:);
t2 = t2(:)';

K = zeros(numel(te), numel(t2));
for j = 1:numel(t2)
    K(:,j) = exp(-te/t2(j));
end
K = [K ones(size(K,1), 1)]; % for noise fitting

% Regularization matrix
R = eye(size(K,2)); 
for ii=1:size(R,1)-1, 
    R(ii,ii+1) = -1; 
end;

KR = [K; sqrt(lambda)*R];

f = lsqnonneg(KR, [y; zeros(size(K,2),1)]);
y_fit = K*f;

b = f(end); % the estimated Rician noise mean

figure(4),
subplot(1,2,1),
plot(te, y, 'b o'), hold on
xlabel('TE (ms)'), ylabel('Signal Intensity'), title('T2 decay curve'),
plot(te, y_fit, 'r*'),
h = line([te(1) te(end)], [b(end) b(end)]); set(h,'color','black');
legend('Measured','Fitted','Base line')

subplot(1,2,2),
semilogx(t2, f(1:end-1))
xlabel('T2 (ms)'), title('T2 spectrum'), grid on
axis([t2(1) t2(end) 0 1.1*max(f(:))])
