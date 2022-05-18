function [ Tu,count ] = AP_Gridless( Y,r,K,max_iter,tol,Tu_init,make_plot )
%Mark Wagner, Dec 31 2019
%extended Alternating Projection based gridless beamforming (e-APG) for
%arbitrary array geometry. Given measurements Y, and number of sources K,
%estimate DOAs. max_iter is maximum number of iterations, tolerance is a
%break condition calling the function to end angle between subspaces of Tu
%and the projection of its extended Vandermonde matrix is low enough.

% M = number of sensors
% L = number of snapshots
% K = number of sources
% Y = MxL measurement matrix
% r = vector of sensor positions
% Tu= MxM hermitian toeplitz matrix of rank K
% Z = LxL Hermitian symmetric matrix
% max_iter = maximum number of iterations
% tol = break condition

%% end initialization
Y           = (Y./norm(Y,'fro'));
r           = r(:).';
[M,L]       = size(Y);                  %number of sensors, number of snapshots

Tu              = Tu_init;
Z               = diag(diag(Y'*Y)/M); 
Bs              = zeros(max_iter,1);        %allocate space to record break conditions
S_prev          = zeros(M+L);
count           = 0;
for i = 1:max_iter
    S               = [Tu,Y;Y',Z];      %construct S
    S               = S./norm(S);
    S               = nearestSPD(S);    %PSD projection
    Tu              = S(1:M,1:M);
    [Tu,~,~]        = ETP(K,r,Tu);      % Extended Toeplitz Projection of matrix Tu

    Bs(i)               = norm(S-S_prev,'fro') / norm(S_prev,'fro');
%     Bs(i)               = norm(S-S_prev,'fro');
    if Bs(i) < tol
        count = i;
        break
    end
    S_prev = S;    
end

if count == 0
    count = max_iter;
end
if make_plot
    figure(3)
    loglog(abs((Bs)))
    grid on
    title('$\| S^{(k)}- S^{(k+1)} \|_F$ vs. Iteration', 'Interpreter','latex')
    xlabel('Iteration')
    ylabel('$\| S^{(k)}- S^{(k+1)} \|_F$', 'Interpreter','latex')
    pause(.01)
end
end

function Ahat = nearestSPD(A)
% nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
% usage: Ahat = nearestSPD(A)
%
% From Higham: "The nearest symmetric positive semidefinite matrix in the
% Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
% where H is the symmetric polar factor of B=(A + A')/2."
%
% http://www.sciencedirect.com/science/article/pii/0024379588902236
%
% arguments: (input)
%  A - square matrix, which will be converted to the nearest Symmetric
%    Positive Definite Matrix.
%
% Arguments: (output)
%  Ahat - The matrix chosen as the nearest SPD matrix to A.
if nargin ~= 1
  error('Exactly one argument must be provided.')
end
% test for a square matrix A
[r,c] = size(A);
if r ~= c
  error('A must be a square matrix.')
elseif (r == 1) && (A <= 0)
  % A was scalar and non-positive, so just return eps
  Ahat = eps;
  return
end
% symmetrize A into B
B = (A + A')/2;
% Compute the symmetric polar factor of B. Call it H.
% Clearly H is itself SPD.
[~,Sigma,V] = svd(B);
H = V*Sigma*V';
% get Ahat in the above formula
Ahat = (B+H)/2;
% ensure symmetry
Ahat = (Ahat + Ahat')/2;
% test that Ahat is in fact PD. if it is not so, then tweak it just a bit.

p = 1;
k = 0;
count = 0;
while p ~= 0 && count < 10000
  [R,p] = chol(Ahat);
  k = k + 1;
  if p ~= 0
    % Ahat failed the chol test. It must have been just a hair off,
    % due to floating point trash, so it is simplest now just to
    % tweak by adding a tiny multiple of an identity matrix.
    mineig = min(eig(Ahat));
    Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
  end
  count = count+1;
end
end

function [ Y, V, D ] = ETP(K, r, X )
%Extended Toeplitz Projection of matrix X

[t_est,c1]          = vandermonde(r, K, X );    %decompose
alpha_est           = exp(1i*2*pi*t_est);       %reconstruct alpha parameters
D                   = (diag(c1));               %signal strengths
V                   = (alpha_est(:).').^r(:);   %reconstruct extended Vandermonde matrix
Y                   = (V*D*V');                 %reconstruct Tu

end

function [ root_locs, c1 ] = vandermonde( r, K, T )
% Mark Wagner, Dec 31, 2019
% GENERALIZED VANDERMONDE DECOMPOSITION: FINDS DOAS FROM COVARIANCE MATRIX FORMED FROM
% UNEVENLY SAMPLED SUM OF SIGNALS IN NOISE 
% INPUTS: 
% r             SENSOR LOCATIONS 
% K             NUMBER OF DOAS PRESENT
% T             COVARIANCE MATRIX
% OUTPUT:
% root_locs     DOAS GIVEN BETWEEN [-.5, .5)
% c1            SOURCE POWERS

%% Make sure inputs are correct size, find # sensors
r               = real(r(:)).';
M               = size(T,1);

%% ULA case
if diag(squareform(pdist(r')),1) == ones(M-1,1)
    root_locs   = rootmusic(T,K,'corr');
    root_locs   = root_locs./(2*pi);
    root_locs   = root_locs(:);
    W_est       = exp(1i*2*pi*root_locs'.*r');%regenerate irregular Vandermonde Matrix 
    c1          = real(diag(pinv(W_est)*T*pinv(W_est)'));
    return
end

%% Find noise subspace
[U,~,~]         = svd(T);
Un              = U(:,K+1:end);
G               = Un*Un';

%% Evaluate MUSIC spectrum
samples                 = 100*(M);
spacing                 = 180/samples;
f                       = -90:spacing:90;
% f                       = [-90-spacing,f,90+spacing];
MUSIC_spectrum          = zeros(length(f),1);
for i = 1:length(f)
    a                   = exp(1i*pi*sind(f(i)).*r');        %steering vector
    MUSIC_spectrum(i)   = -abs((a'*G*a));               %negative null spectrum 
end

%% find peaks
[pks,inds]              = findpeaks(MUSIC_spectrum,(1:length(MUSIC_spectrum))); %find the peaks
[~,id]                  = sort(pks,'descend');

%% refine root estimates to high precision
root_locs               = f(inds(id(1:min([length(id),K]))));
root_locs_refined       = zeros(K,1);
% fun                     = @(f) abs(( (exp(1i*pi*sind(f).*r(:).'))*G'*exp(1i*pi*sind(-f).*r(:)) ));     %generalized null spectrum 
fun                     = @(f) abs(diag((exp(1i*pi*sind(f).*r')')*G*exp(1i*pi*sind(f).*r')));     %generalized null spectrum 
% fun                     = @(f) null_spec_polynomial(G, r, exp(1i*pi*sind(f)) );

options                 = optimset('Display','none','TolX',1e-25);
for i = 1:length(root_locs)                                                             %for each approximate root
%     root_locs_refined(i)= fminsearch(fun,root_locs(i)-(abs(f(1)-f(2))/2),options);      %find local minima about generalized null spectrum
    root_locs_refined(i)= fminbnd(fun,root_locs(i)-(abs(f(1)-f(2))/2),root_locs(i)+(abs(f(1)-f(2))/2),options);      %find local minima about generalized null spectrum
end      
rl_rad                  = root_locs_refined/180*pi;
root_locs               = sort(sin(rl_rad)/2);

%% Reconstruct

W_est                   = exp(1i*2*pi*root_locs(:).'.*r(:));%regenerate irregular Vandermonde Matrix
W_inv                   = pinv(W_est);
c1                      = real(diag(W_inv*T*W_inv'));

end

