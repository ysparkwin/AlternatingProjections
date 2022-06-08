function [ Tu,count,Tu1st,Tu2nd,Tu3rd ] = APCOVv1p00(Y,q,Nalg,max_iter,tol,Tu_init,Z_init)

% number of sensors, number of snapshots
[My,Ly] = size(Y);
% sample covariance matrix (SCM)
Rhat = Y*Y' / Ly;
% Normalization for making the initialization have a consistent energy
Rhat = Rhat/norm(Rhat,'fro')*My*My;

% Initialization
Tu_init = Tu_init/norm(Tu_init,'fro')/1;
Z = Z_init;
Z = Z/norm(Z,'fro')/1;

count           = 0;
Tu              = Tu_init;
Bs              = zeros(max_iter,1);        %allocate space to record break conditions
S_prev          = zeros(My+My);
for i = 1:max_iter
    
    S = [Tu Rhat; Rhat Z];
    S = S./norm(S,'fro')*My*My;   % Normalization for energy conservation
    S = PSD(S);             % PSD projection
    if i==1, Tu1st = Tu; end
    if i==2, Tu2nd = Tu; end
    if i==3, Tu3rd = Tu; end
    
    Tu = S(1:My,1:My);
    [Tu,~,~] = ETP(Nalg,q,Tu); % Extended Toeplitz Projection of matrix Tu

    Z                   = S(My+1:end,My+1:end);
    if i==1, Tu1st = Tu; end
    if i==2, Tu2nd = Tu; end
    if i==3, Tu3rd = Tu; end
    Bs(i)               = norm(S-S_prev,'fro') / norm(S_prev,'fro');
%     Bs(i)               = norm(S-S_prev,'fro');
%     Bs(i)               = norm(S(1:My,1:My)-S_prev(1:My,1:My),'fro');
    
    if Bs(i) < tol
        count =         i;
        if exist('Tu2nd','var') == 0, Tu2nd = Tu; end
        if exist('Tu3rd','var') == 0, Tu3rd = Tu; end
        break
    end
    S_prev = S;
    
end
if count == 0
    count = max_iter;
end
end

function [ Y ] = PSD( X )
%Postive semi-definite (PSD) projection of a matrix
[VZ,DZ]             = eig((X+X')/2);
dZ                  = real(diag(DZ));                   %eigenvalues must be real
idx                 = dZ>0;
Y                   = VZ(:,idx)*diag(dZ(idx))*VZ(:,idx)';%reconstruct
Y                   = (Y+Y')/2;
end

function [ Y, V, D ] = ETP(K, r, X )
%Extended Toeplitz Projection of matrix X

[t_est,c1]          = vandermonde(r, K, X );    %decompose
alpha_est           = exp(1i*2*pi*t_est);       %reconstruct alpha parameters
D                   = (diag(c1));               %signal strengths
V                   = (alpha_est(:).').^r(:);   %reconstruct extended Vandermonde matrix
Y                   = (V*D*V');                 %reconstruct Tu

% if sum(isnan(Y(:)))+ sum(isinf(Y(:))) > 0
%     Y = zeros(size(X));
% end

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

