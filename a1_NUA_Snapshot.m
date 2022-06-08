% Version 1.0: (05/09/2022)
% written by Y. Park
% System framework as defined in Noiselab DOA estimation

% Version 2.0: (06/08/2022)
% written by Y. Park
% AP Covariance (JASA 2022)

% Mark Wagner, Yongsung Park, & Peter Gerstoft
% MPL/SIO/UCSD
% wagnermark1992@gmail.com / yongsungpark@ucsd.edu / gerstoft@ucsd.edu
% noiselab.ucsd.edu

% Citation
% M. Wagner, Y. Park, and P. Gerstoft, “Gridless DOA estimation and root-MUSIC for non-uniform linear arrays,” IEEE Trans. Signal Process. 69, 2144–2157 (2021).
% M. Wagner, P. Gerstoft, and Y. Park, “Gridless DOA estimation via alternating projections,” in Proc. IEEE ICASSP (2019), pp. 4215–4219.
% Y. Park and P. Gerstoft, “Alternating pro jections gridless covariance-based estimation for DOA,” in Proc. IEEE ICASSP (2021), pp. 4385–4389.
% Y. Park and P. Gerstoft, “Gridless sparse covariance-based beamforming via alternating projections including co-prime arrays,” J. Acoust. Soc. Am. 151(6),3828-3837 (2022), pp. 4385–4389.

%%
clear; clc; %close all;

addpath([cd,'/_common'])
errCut = 10; % Maximum RMSE cut-off.

dbstop if error;

% Case: NUA       Multi- snapshot (for NUA, see Line 47-49)
SNR = 20;
% anglesTrue = [-65; -2; 3]; % +- .5 Line 73

xAxes = linspace(0.4771,2,6);
xAxes = round(10.^xAxes); % Number of Snapshots list
for xaxis = 1:length(xAxes)
% Number of Monte-Carlo simulations
Nsim = 100;
for nsim=1:Nsim

rngN = (xaxis-1)*Nsim + nsim; rng(rngN);
disp(' ')
disp(['Number of Snapshots',num2str(xAxes(xaxis)),'#Sim : ',num2str(nsim)])

% Environment parameters
c = 1500;       % speed of sound
f = 200;        % frequency
lambda = c/f;   % wavelength

% ULA-horizontal array configuration
Nsensor = 20;               % number of sensors
d = 1/2*lambda;             % intersensor spacing
q = (0:1:(Nsensor-1))';     % sensor numbering

% sensor spacing perturbation for NUA
qPert = rand(size(q))-.5; qPert(1) = 0; qPert(end) = 0;
q = q+qPert;

xq = q*d;   % sensor locations

% signal generation parameters
% SNR = xAxes(xaxis);

% total number of snapshots
Nsnapshot = xAxes(xaxis);

% range of angle space
thetalim = [-90 90];

% Angular search grid
theta_separation = .25;
theta = (thetalim(1):theta_separation:thetalim(2))';
Ntheta = length(theta);

% Design/steering matrix (Sensing matrix)
sin_theta = sind(theta);
sensingMatrix = exp(-1i*2*pi/lambda*xq*sin_theta.')/sqrt(Nsensor);

% Generate received signal
anglesTrue = [-65; -2; 3];
anglesTrue = anglesTrue + rand(size(anglesTrue)) - 0.5;

disp(['True DOAs : ',num2str(anglesTrue.')])
% source_amp = [  7;   7;  7;  4;  4; 13];
anglesTracks = repmat(anglesTrue,[1,Nsnapshot]);
sinAnglesTracks = sind(anglesTracks);
Nsource = numel(anglesTrue);

receivedSignal = zeros(Nsensor,Nsnapshot);
for snapshot = 1:Nsnapshot
    % Source amplitude
%     source_amp(:,snapshot) = 6*rand(size(anglesTrue)) + 4;
    source_amp(:,snapshot) = 10*ones(size(anglesTrue));
%     source_amp(:,snapshot) = [10; 7; 4; 7];
    Xsource = source_amp(:,snapshot).*exp(1i*2*pi*rand(Nsource,1));    % random phase
    
    % Represenation matrix (steering matrix)
    transmitMatrix = exp( -1i*2*pi/lambda*xq*sinAnglesTracks(:,snapshot).' );
    
    % Received signal without noise
    receivedSignal(:,snapshot) = sum(transmitMatrix*diag(Xsource),2);
    
    % add noise to the signals
    rnl = 10^(-SNR/20)*norm(Xsource);
    nwhite = complex(randn(Nsensor,1),randn(Nsensor,1))/sqrt(2*Nsensor);
    e = nwhite * rnl;	% error vector
    receivedSignal(:,snapshot) = receivedSignal(:,snapshot) + e;
    
    % for CRB
    crnl(snapshot) = rnl;
    cX(:,snapshot) = Xsource;
end


%% CRB
Acrb = exp( -1i*2*pi/lambda*xq*sinAnglesTracks(:,snapshot).' );
AcrbD = (-1i*2*pi/lambda*xq*cosd(anglesTracks(:,snapshot)).') ...
    .* exp( -1i*2*pi/lambda*xq*sinAnglesTracks(:,snapshot).' );

Xs = cX./(crnl/crnl(1));
Pn = power(crnl(1),2)/Nsensor;
Phat = diag(diag(Xs*Xs'/Nsnapshot));

H = AcrbD'...
    *(eye(Nsensor) - Acrb/(Acrb'*Acrb)*Acrb')...
    *AcrbD;

% det. CRB
CRB = real(H .* (Phat.'));
CRB = eye(size(Xs,1)) / CRB * (Pn / Nsnapshot / 2);
outputsCRBd(xaxis,nsim) = mean(diag(CRB));

%% Conventional beamforming (CBF)
Ryy = receivedSignal*receivedSignal' / Nsnapshot;
Pcbf = zeros(numel(theta),1);
for ii=1:length(theta)
    Pcbf(ii) = sensingMatrix(:,ii)'*Ryy*sensingMatrix(:,ii)/(sensingMatrix(:,ii)'*sensingMatrix(:,ii));
end
% Pcbf = sensingMatrix' * receivedSignal;
% plot(theta,mean(Pcbf.*conj(Pcbf),2)/max(mean(Pcbf.*conj(Pcbf),2)),'k:','linewidth',1,'displayname','CBF')
% plot(theta,abs(Pcbf)/max(abs(Pcbf)),'k:','linewidth',1.5,'displayname','CBF')

[~, Ilocs] = findpeaks(abs(Pcbf),'SORTSTR','descend','Npeaks', Nsource);
% DoA_error = errorDOA(theta(Ilocs),anglesTrue);
DoA_error = errorDOAcutoff(theta(Ilocs),anglesTrue,errCut);
disp(['RMSE CBF: ',num2str(sqrt(mean(power(DoA_error,2))))])

if nsim==1 && xaxis==1, outputsCBF = []; end
outputCBF = struct('theta',theta(Ilocs),'error',DoA_error);
outputsCBF = [outputsCBF; outputCBF];


%% root-MUSIC, "NUA", irregular root-MUSIC
[t_est,~]   = wagner_decomp( q, Nsource, Ryy );
t_est = -t_est*lambda/d;
DoA_est_deg = asin(t_est)/pi*180;
DoA_error = errorDOAcutoff(DoA_est_deg,anglesTrue,errCut);
disp(['RMSE irr. root-MUSIC: ',num2str(sqrt(mean(power(DoA_error,2))))])

% %for ULA
% spRmusic = rmusic_1d(Ryy, Nsource, 2*pi*d/lambda);
% DoA_error = errorDOAcutoff(-rad2deg(spRmusic.x_est),anglesTrue,errCut);
% disp(['RMSE root-MUSIC: ',num2str(sqrt(mean(power(DoA_error,2))))])

if nsim==1 && xaxis==1, outputsrMUSIC = []; end
% outputrMUSIC = struct('theta',-rad2deg(spRmusic.x_est),'error',DoA_error);
outputrMUSIC = struct('theta',DoA_est_deg,'error',DoA_error);
outputsrMUSIC = [outputsrMUSIC; outputrMUSIC];


%% AP-ULA, "NUA"

% MAX_IT = 1000;
% NeigAP = Nsource;
% tol    = 1e-3;
% 
% [T,iAPM]    = AP_Gridless( receivedSignal,q,NeigAP,MAX_IT,tol,zeros(Nsensor),0 );
% [t_est,~]   = wagner_decomp( q, NeigAP, T );     %decompose
% t_est = -t_est*lambda/d;
% DoA_est_deg = asin(t_est)/pi*180;
% DoA_error = errorDOAcutoff(DoA_est_deg,anglesTrue,errCut);
% disp(['RMSE AP-Snapshot ULA: ',num2str(sqrt(mean(power(DoA_error,2))))])
% 
% if nsim==1 && xaxis==1, outputsAPula = []; end
% outputAPula = struct('theta',DoA_est_deg,'error',DoA_error);
% outputsAPula = [outputsAPula; outputAPula];


%% AP-Gridless for NUA

MAX_IT = 1000;
NeigAP = Nsource;
tol    = 1e-3;

% [T,iAPM]    = AP_Gridless( receivedSignal,q,NeigAP,MAX_IT,tol,zeros(Nsensor),0 );  % for AP-ULA
[T,iAPM]    = AP_Gridless( receivedSignal,q+.1,NeigAP,MAX_IT,tol,zeros(Nsensor),0 );    % for AP-Gridless
[t_est,~]   = wagner_decomp( q, NeigAP, T );     %decompose
t_est = -t_est*lambda/d;
DoA_est_deg = asin(t_est)/pi*180;
DoA_error = errorDOAcutoff(DoA_est_deg,anglesTrue,errCut);
disp(['RMSE AP-Snapshot: ',num2str(sqrt(mean(power(DoA_error,2))))])

if nsim==1 && xaxis==1, outputsAPsnapshot = []; end
outputAPsnapshot = struct('theta',DoA_est_deg,'error',DoA_error);
outputsAPsnapshot = [outputsAPsnapshot; outputAPsnapshot];


%% AP-Covariance

max_iter = 1000;
Nalg    = Nsource;
tol     = 1e-4;

if exist('Tu_init','var') == 0
    Tu_init = rand(Nsensor) + 1i*rand(Nsensor);
    Z_init = rand(Nsensor) + 1i*rand(Nsensor);
end
[ Tu,iAP,~,~,~ ]  = APCOVv1p00(receivedSignal,q,Nalg,max_iter,tol,Tu_init,Z_init);

[t_est,~]   = wagner_decomp( q, Nalg, Tu );     %decompose
t_est = -t_est*lambda/d;
while(1)
    t_est(t_est>1) = t_est(t_est>1)  - 2;
    t_est(t_est<-1)= t_est(t_est<-1) + 2;
    if sum(t_est>1 | t_est<-1) == 0, break; end
end
DoA_est_deg = asin(t_est)/pi*180;
clear Tu_init Z_init
DoA_error = errorDOAcutoff(DoA_est_deg,anglesTrue,errCut);
disp(['RMSE AP-Covariance: ',num2str(sqrt(mean(power(DoA_error,2))))])

if nsim==1 && xaxis==1, outputsAPcov = []; end
outputAPcov = struct('theta',DoA_est_deg,'error',DoA_error);
outputsAPcov = [outputsAPcov; outputAPcov];


%%  SBL
options = SBLSet();
options.convergence.error = 10^(-3);
options.Nsource = ceil(Nsensor/2);
options.gamma_range=10^-20;

[gamma, reportSBL] = SBL_v4( sensingMatrix, receivedSignal, options );
[~,peak_SBL] = findpeaks(gamma,'SORTSTR','descend','Npeaks', Nsource);
DoA_error = errorDOAcutoff(theta(peak_SBL),anglesTrue,errCut);
disp(['RMSE SBL: ',num2str(sqrt(mean(power(DoA_error,2))))])

if nsim==1 && xaxis==1, outputsSBL = []; end
outputSBL = struct('theta',theta(peak_SBL),'error',DoA_error);
outputsSBL = [outputsSBL; outputSBL];


end

end



%% Figure
for ind=1:length(xAxes)
totETcbf = [];
totETrmu=[];
% totETAPula=[];
totETAPsnapshot=[];
totETAPcov=[];
totETsbl = [];
for index=1:Nsim
    totETcbf = [totETcbf;outputsCBF((ind-1)*Nsim+index).error];
    totETrmu = [totETrmu;outputsrMUSIC((ind-1)*Nsim+index).error];
%     totETAPula = [totETAPula;outputsAPula((ind-1)*Nsim+index).error];
    totETAPsnapshot = [totETAPsnapshot;outputsAPsnapshot((ind-1)*Nsim+index).error];
    totETAPcov = [totETAPcov;outputsAPcov((ind-1)*Nsim+index).error];
    totETsbl = [totETsbl;outputsSBL((ind-1)*Nsim+index).error];
end

Nout = 0; % Portion of Outliers, (ignore)
totETcbf = sort(abs(totETcbf));
totETrmu = sort(abs(totETrmu));
% totETAPula = sort(abs(totETAPula));
totETAPsnapshot = sort(abs(totETAPsnapshot));
totETAPcov = sort(abs(totETAPcov));
totETsbl = sort(abs(totETsbl));

ecbf(ind) = sqrt(mean(power(totETcbf(1:length(totETcbf)-floor(length(totETcbf)*Nout)),2)));
ermusic(ind) = sqrt(mean(power(totETrmu(1:length(totETrmu)-floor(length(totETrmu)*Nout)),2)));
% eapula(ind) = sqrt(mean(power(totETAPula(1:length(totETAPula)-floor(length(totETAPula)*Nout)),2)));
eapsnapshot(ind) = sqrt(mean(power(totETAPsnapshot(1:length(totETAPsnapshot)-floor(length(totETAPsnapshot)*Nout)),2)));
eapcov(ind) = sqrt(mean(power(totETAPcov(1:length(totETAPcov)-floor(length(totETAPcov)*Nout)),2)));
esbl(ind) = sqrt(mean(power(totETsbl(1:length(totETsbl)-floor(length(totETsbl)*Nout)),2)));
end

figure; set(gcf,'position',[750,200,700,600]);
hold on;
h1=plot(xAxes,sqrt(mean(outputsCRBd,2)*180/pi*180/pi),'k','linewidth',1.0,'markersize',10,'displayname','det. CRB');
figH = h1;
h3=plot(xAxes,ecbf,'k:','linewidth',1.5,'markersize',10,'displayname','CBF');
figH = [figH,h3];
pcolor = lines;
h5=plot(xAxes,esbl,'m--','linewidth',1.2,'markersize',10,'displayname',...
    'SBL','color',pcolor(2,:)); figH = [figH,h5];
h8=plot(xAxes,ermusic,'b-.','linewidth',2,'markersize',10,'displayname',...
    'irr. Root-MUSIC'); figH = [figH,h8];
% h9=plot(xAxes,eapula,'c-.','linewidth',1.8,'markersize',12,'displayname',...
%     'AP-ULA'); figH = [figH,h9];
h91=plot(xAxes,eapsnapshot,'linewidth',1.8,'markersize',12,'color',pcolor(7,:),'displayname',...
    'AP-Snapshot'); figH = [figH,h91];
h92=plot(xAxes,eapcov,'linewidth',1.8,'markersize',12,'color','r','displayname',...
    'AP-Covariance'); figH = [figH,h92];
hold off;

% xlabel('SNR~[dB]','interpreter','latex')
xlabel('Number of snapshots','interpreter','latex')
ylabel('RMSE~[$^\circ$]','interpreter','latex')
legend(fliplr(figH),'location','northeast','interpreter','latex')

set(gca,'fontsize',18,'yscale','log','xscale','log')
box on; grid on;

axis([min(xAxes) max(xAxes) min(sqrt(mean(outputsCRBd,2)*180/pi*180/pi)) 1])
% axis([0 20 0 10])

%%
rmpath([cd,'/_common'])
