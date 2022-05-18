function [ root_locs, c1 ] = wagner_decomp( r, K, T )
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
samples                 = 100*(M);%15*(max(r)-min(r));%10*(M);
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