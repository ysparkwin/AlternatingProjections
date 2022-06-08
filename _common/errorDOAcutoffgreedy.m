function [err] = errorDOAcutoff(theEst,theTrue,cutoff)
err = zeros(size(theTrue));

if size(theEst,1)<size(theEst,2), theEst = theEst.'; end
if size(theTrue,1)<size(theTrue,2), theTrue = theTrue.'; end

theE = theEst;
if numel(theEst)<numel(theTrue)
    theE = [theE;-inf*ones(ceil(numel(theTrue)/numel(theEst))-1,1)];
end
theT = theTrue;

Ap = (repmat(theT',numel(theE),1)-repmat(theE,1,numel(theT)));
for i=1:numel(theTrue)
    [ind1,ind2] = find(abs(Ap)==min(min(abs(Ap))));
%     err(i) = Ap(ind1(1),ind2(1));
    err(i) = min(abs(Ap(ind1(1),ind2(1))),cutoff);
    
    Ap(ind1(1),:) = []; Ap(:,ind2(1)) = [];
end