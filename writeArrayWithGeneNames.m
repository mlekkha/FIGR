%
% This function is a bit deprecated.
% Ideally we would supply "saveStdGeneExprFiles.m"
%

function writeArrayWithGeneNames (filename, a, geneNames)
delim = {'\t', '\n', '\n', '\n', '\n', '\n', '\n', '\n', '\n', 'n', 'n', 'n'};
fmtstr = '%.2f';
dmax = ndims(a);
cmaxd = size(a);
fid = fopen (filename, 'w');
for i=1:numel(geneNames)
    fprintf (fid, '%s', string(geneNames(i)));
    fprintf (fid, '\t');
end
fprintf (fid, '\n');
idx = ones (cmaxd);                 % start all indices at 1
for i=1:numel(a)
    fprintf (fid, fmtstr, a(i));      % print array element
    for d=1:dmax
        fprintf (fid, delim{d});    % print delimiter
        idx(d) = idx(d)+1;          % increment index
        if (idx(d) <= cmaxd(d)); break; end;
        idx(d) = 1;                 % restart index at 1
    end                             % consider index at next level
end
fclose (fid);
end