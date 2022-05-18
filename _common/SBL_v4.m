function [ gamma, report] = SBL_v4( A , Y, options )
%
% function [ gamma , report ] = SBL_v3p1( A , Y , Nfreq, Nsource,flag )
% The idea behind SBL is to find a diagonal replica 'covariance' Gamma.
% Minimizing (YY^T / AGA^T + penality) should lead to the correct
% replica selection (up to a bogus scale factor/amplitude).
%
% Inputs
%
% A - Multiple frequency augmented dictionary <f , n , m>
%     f: number of frequencies
%     n: number of sensors
%     m: number of replicas
%   Note: if f==1, A = < n , m >
%
%
% Y - Multiple snapshot multiple frequency observations <f , n , L>
%     f: number of frequencies
%     n: number of sensors
%     L: number of snapshots
%   Note: if f==1, Y = < n , L >
%
% options - see SBLset.m 
%
%
% Outputs
%
% gamma <m , 1> - vector containing source power
%                 1: surfaces found by minimum error norm
%
% report - various report options
%
%--------------------------------------------------------------------------
% Version 1.0:
% Code originally written by P. Gerstoft.
%
% Version 2.23
% Edited to include multiple frequency support: 5/16/16
%
% Version 3.1
% Different convergance norm and code update
% A and Y have now one more dimensions
% Posterior unbiased mean
% handles single snapshot
%
% Version 3.32
% Robust version of SBL with Phi parameter in options structure
%
% TODO:
% posterior covariance
%
% Santosh Nannuru & Kay L Gemba
% MPL/SIO/UCSD gemba@ucsd.edu & snannuru@ucsd.edu

% number of frequencies
Nfreq     = size(A,3);
% single frequency single snapshot
%%
options.SBL_v = '4.00';

%% slicing
% 
% if options.tic == 1
%     tic
% end

%% Initialize variables
Nsource = options.Nsource;
% number of sensors
Nsensor   = size(A,1);
% number of dictionary entries
Ntheta    = size(A,2);
% number of snapshots in the data covariance
Nsnapshot = size(Y,2);
% noise power initialization
sigc      = ones(Nfreq,1) * options.noisepower.guess;
%sigc      = options.noisepower.guess;
% posterior
%x_post    = zeros( Ntheta, Nsnapshot,Nfreq);
% minimum (global) gamma
gmin_global = realmax;
% reduce broacast
%phi         = options.phi;
% space allocation
errornorm   = zeros(options.convergence.maxiter,1);

% initialize with CBF output (assume single frequency)
%gamma        = Bartlett_processor(squeeze(A), squeeze(Y));
 gamma        = zeros(Ntheta, 1);
 for iF = 1:Nfreq
    Af = squeeze(A(:,:,iF));
    Yf = squeeze(Y(:,:,iF));
    gamma = gamma+sum(abs(Af' * Yf).^2, 2) / Nsnapshot;
end

%gamma        = ones(Ntheta, 1);
gamma_num    = zeros(Ntheta,Nfreq);
gamma_denum  = zeros(Ntheta,Nfreq);

% Sample Covariance Matrix
SCM = zeros(  Nsensor , Nsensor, Nfreq);
for iF = 1 : Nfreq
    SCM(:,:,iF) = squeeze(Y(:,:,iF)) * squeeze(Y(:,:,iF))' / Nsnapshot;
end

%% Main Loop

display(['SBL version ', options.SBL_v ,' initialized.']);

% override initialization
%gamma = ( sum ( abs ( ( squeeze(A)' * squeeze(Y) ).^2 ),2 )) / Nsnapshot;

for j1 = 1 : options.convergence.maxiter
        
    % for error analysis
    gammaOld = gamma;
    Itheta= find(gamma>max(gamma)*options.gamma_range);
    gammaSmall=gamma(Itheta);
    %% gamma update
    % this is a sum over frequency and can be done by multiple processors
    % --> multi-frequency SBL should be almost as fast as single frequency
    % if num(proc) = Nfreq
    
    for iF = 1 : Nfreq
        Af = squeeze (A(:,Itheta,iF));
        Yf = squeeze(Y(:,:,iF));
        % SigmaY inverse
        %SigmaYinv = eye(Nsensor) / (sigc(i_f) * eye(Nsensor) + ...
        %   Af * (repmat(gamma, [1 Nsensor] ) .* Af') + phi * sum(gamma) * eye(Nsensor));
        SigmaYinv = eye(Nsensor) / (sigc(iF) * eye(Nsensor) + ...
           Af * bsxfun(@times, gammaSmall, Af') );
 %         Af * bsxfun(@times, gammaSmall, Af') + phi * sum(gamma) * eye(Nsensor));

%         SigmaYinvY = SigmaYinv * Yf;
        
%         % Sum over snapshots and normalize
%         gamma_num(Itheta,iF)   = ( sum ( abs ( ( Af' * SigmaYinv * Yf ).^2 ),2 ) + phi * abs(SigmaYinvY(:)' * SigmaYinvY(:)) ) / Nsnapshot;
% 
%         % positive def quantity, abs takes care of roundoff errors        
%         gamma_denum(Itheta,iF) = abs( sum  ( (Af' * SigmaYinv).' .* Af, 1 ) ) + phi * abs(sum(diag(SigmaYinv))) ;
        % Sum over snapshots and normalize
        gamma_num(Itheta,iF)   = ( sum ( abs ( ( Af' * SigmaYinv * Yf ).^2 ),2 ) )  / Nsnapshot;

        % positive def quantity, abs takes care of roundoff errors        
        gamma_denum(Itheta,iF) = abs( sum  ( (Af' * SigmaYinv).' .* Af, 1 ) )  ;

    end
    
    % Fixed point Eq. update
    gamma(Itheta)  = gamma(Itheta)   .* ((sum( gamma_num(Itheta,:)    ,2 ) ./...
        sum( gamma_denum(Itheta,:)  ,2 ) ).^(1/options.fixedpoint) ) ;
    %% sigma update
    
%     locate same peaks for all frequencies
%     [ ~ , Ilocs] = findpeaks(gamma,'SORTSTR','descend','NPEAKS',Nsource);
    [~, Ilocs] = SBLpeaks_1D(gamma, Nsource);

    for iF = 1 : Nfreq
        Am        = squeeze(A(:,Ilocs,iF));        % only active replicas
        % noise estimate
        sigc(iF) = real(trace( (eye(Nsensor)-Am*pinv(Am)) * squeeze(SCM(:,:,iF)) ) / ( Nsensor - Nsource ) );
    end
    
    %plot(gamma); grid on; pause(0.005);
    
    %% Convergence  checks convergance and displays status reports
    % convergence indicator
    errornorm(j1) = norm ( gamma - gammaOld, 1 ) / norm ( gamma, 1 );
        
    % global min error
    if j1 > options.convergence.min_iteration  &&  errornorm(j1) < gmin_global
        gmin_global  = errornorm(j1);
        gamma_min    = gamma;
        iteration_L1 = j1;
    end
    
    % inline convergence code
    if j1 > options.convergence.min_iteration && ( errornorm(j1) < options.convergence.error  || iteration_L1 + options.convergence.delay <= j1)
        if options.flag == 1
            display(['Solution converged. Iteration: ',num2str(sprintf('%.4u',j1)),'. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.'])
        end
        break; % goodbye     
    elseif j1 == options.convergence.maxiter % not convereged
        if options.flag == 1
            warning(['Solution not converged. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.'])
        end
        % status report
    elseif j1 ~= options.convergence.maxiter  && options.flag == 1 && mod(j1,options.status_report) == 0 % Iteration reporting
        display(['Iteration: ',num2str(sprintf('%.4u',j1)),' Dic size: ',num2str(length(Itheta)),'. Error: ',num2str(sprintf('%1.2e' , errornorm(j1) )),'.' ])
    end
    
end

%% Posterior distribution for polarity
% for iF = 1 : Nfreq
%     Af = squeeze (A(:,:,iF));
%     x_post(:,:,iF) = repmat(gamma, [1 Nsnapshot] ) .* (Af' / (sigc(iF) * eye(Nsensor) + ...
%         Af * (repmat(gamma, [1 Nsensor] ) .* Af') + phi * sum(gamma) * eye(Nsensor)) * squeeze(Y(:,:,iF)));
% end

%% function return
gamma = gamma_min;

%% Report section
% vectors containing errors
report.results.error    = errornorm;
% Error when minimum was obtained
report.results.iteration_L1 = iteration_L1;
% General info
report.results.final_iteration.iteration = j1;
report.results.final_iteration.noisepower = sigc;

% debug output parameters (assuming single frequency)
%report.SigmaYinv = SigmaYinv;
%report.SCM = squeeze(SCM);
% 
% if options.tic == 1
%     report.results.toc = toc;
% else
%     report.results.toc = 0;
% end
% data
report.results.final_iteration.gamma  = gamma  ;
%report.results.final_iteration.x_post = x_post ;
report.options = options;
end

function [pks, locs] = SBLpeaks_1D(gamma, Nsources)
%
% [pks, locs] = SBLpeaks_1D(gamma, Nsources)
%
% fast alternative for findpeaks in 1D case
%

% output variables
pks = zeros(Nsources,1);
locs = zeros(Nsources,1);

% zero padding on the boundary
gamma_new = zeros(length(gamma)+2,1);
gamma_new(2:end-1) = gamma;

[~, Ilocs]= sort(gamma,'descend');

% current number of peaks found
npeaks = 0;

for ii = 1:length(Ilocs)
    
    % local patch area surrounding the current array entry i.e. (r,c)
    local_patch = gamma_new(Ilocs(ii):Ilocs(ii)+2);
    
    % zero the center
    local_patch(2) = 0;
    
    if sum(sum(gamma(Ilocs(ii)) > local_patch)) == 3
        npeaks = npeaks + 1;
        
        pks(npeaks) = gamma(Ilocs(ii));
        locs(npeaks) = Ilocs(ii);
        
        % if found sufficient peaks, break
        if npeaks == Nsources
            break;
        end
    end
    
end

% if Nsources not found
if npeaks ~= Nsources
    pks(npeaks+1:Nsources) = [];
    locs(npeaks+1:Nsources) = [];
end

end