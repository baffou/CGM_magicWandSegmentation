% function that extracts the file name 'filename.ext'
% of a string looking like 'path/filename.ext'

function fileName=extractFileName(path)
poss=find(path=='/');
if isempty(poss)
    poss=0;
end
pos=poss(end);
fileName=path(pos+1:end);
