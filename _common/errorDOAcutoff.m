function [err] = errorDOAcutoff(theEst,theTrue,cutoff)

if numel(theTrue) > 10
    disp('Too many true elements for this code.')
    return;
end

if size(theEst,1)<size(theEst,2), theEst = theEst.'; end
if size(theTrue,1)<size(theTrue,2), theTrue = theTrue.'; end

if numel(theEst)<numel(theTrue)
    theEst = [theEst;-inf*ones(ceil(numel(theTrue)/numel(theEst))-1,1)];
end

caseListTmp = perms(1:numel(theEst));
caseList    = unique(caseListTmp(:,1:numel(theTrue)),'rows');

errTmp = theEst(caseList) - repmat(theTrue.',[size(caseList,1),1]);
errTmp(abs(errTmp)>cutoff) = cutoff;

[~, Loc] = min(sum(power(errTmp,2),2));
err = abs(errTmp(Loc(1),:).');

end