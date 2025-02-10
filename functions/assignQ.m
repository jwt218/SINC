function [AQ] = assignQ(CX,IX,varargin)
% assignQ - Matches sample names from file containing age information.
% Inputs:
%   CX       - Output from combineQ [struct]
%   IX       - Output from infoQ [struct]
%              Defaults will be used for SampleMatchLen, Compounds, AdditionalComp, Mode,
%              and AgeFile if not specified.
%   Plot     - Option to generate output figures into a directory ('yes' or 'no')
% Outputs: {AQ] is a struct with fields AQ.
%
% Example Usage:
%   AQ = assignQ(CX, IX, 'Plot', 'yes');

defPlot = 'yes';

expPlot = {'yes','no'};

p = inputParser;
validCX = @(x) isstruct(x);
validIX = @(x) isstruct(x);
validPlot = @(x) any(validatestring(x,expPlot));

addRequired(p,'CX',validCX)
addRequired(p,'IX',validIX)
addParameter(p,'Plot',defPlot,validPlot)

parse(p,CX,IX,varargin{:})

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

CX = p.Results.CX;
IX = p.Results.IX;
Plot = char(p.Results.Plot);
Mode = IX.Mode;

AQ = struct();

tid = IX.SampleMatchLen;
cmp = [IX.Compounds;IX.AdditionalComp];
if isempty(IX.AgeFile)
    uaf = char(unique(table2array(CX.Sample(:,'Identifier1'))));
    uaf2 = unique(string(uaf(:,1:tid)));
    a = table(uaf2, (1:length(uaf2))');
else
    af = IX.AgeFile;
    a = readtable(af);
end

M = CX.Sample;

%if isempty(af); error("Age file must be provided to determine weighted values"); end

dv = M.CorrectedDeltaValue;
comp = M.Component;
lab = M.Identifier;
pks = M.PeakArea;

alab = string(table2array(a(:,1)));
age = table2array(a(:,2));
lalb = char(lab);
ctar = lalb(:,1:tid);
stl = string(ctar);
[ctun,~] = unique(stl,'stable');

ai = ismember(alab,ctun);
tlab = alab(ai);
tage = age(ai);

kage = zeros([length(stl) 1]);
for i = 1:length(tlab)
    ctk = contains(stl,tlab(i));
    kage(ctk) = tage(i);
end

dfa = [kage dv comp pks];
dfs = sortrows(dfa,3); dfs3 = dfs(:,3); udfs3 = unique(dfs3);
DFC = []; pkk = zeros([length(udfs3) 1]);
for i = 1:length(udfs3)
    DFC(i).X = dfs((udfs3(i) == dfs3),:);
    pkk(i) = mean(DFC(i).X(:,4));
end

sz = [length(ctun) length(cmp)];
dxi = zeros(sz); sxi = zeros(sz); pxi = zeros(sz);
for j = 1:length(ctun)

    ki = string(ctar) == string(ctun(j));
    cmpki = comp(ki); dki = dv(ki); pki = pks(ki); Ca = [];
    for i = 1:numel(cmp)
        cd = find(cmpki == cmp(i));
        dm = mean(dki(cd));
        cmm = mean(cmpki(cd));
        smm = std(dki(cd),'omitnan');
        pmm = sum(pki(cd));
        if cd < 2; smm = NaN; end
        Ca = [Ca; dm cmm smm pmm];
    end
    dx = Ca(:,1); dxi(j,:) = dx'; sx = Ca(:,3); sxi(j,:) = sx';
    px = Ca(:,4); pxi(j,:) = px';

end

dz = [length(ai) length(cmp)];
vi = string(cmp);
jri = ["ID";"Age";vi];
tia = find(ai == 1);
dri = NaN(dz); sri = NaN(dz); pri = NaN(dz);

for i = 1:length(tia); dri(tia(i),:) = dxi(i,:); end
dci = [repmat(9999, dz(1), 1) repmat(9999, dz(1), 1) dri];
D = array2table(dci);
D.Properties.VariableNames = jri;
D.ID = alab; D.Age = age;

for i = 1:length(tia); sri(tia(i),:) = sxi(i,:); end
sci = [repmat(9999, dz(1), 1) repmat(9999, dz(1), 1) sri];
S = array2table(sci);
S.Properties.VariableNames = jri;
S.ID = alab; S.Age = age;

for i = 1:length(tia);pri(tia(i),:) = pxi(i,:); end
pci = [repmat(9999, dz(1), 1) repmat(9999, dz(1), 1) pri];
P = array2table(pci);
P.Properties.VariableNames = jri;
P.ID = alab; P.Age = age;

AQ.MeanDelta = D;
AQ.SDDelta = S;
AQ.PeakArea = P;
AQ.Sample = CX.Sample;
AQ.Standard = CX.Standard;
AQ.Age = age;
AQ.CorrectionStds = CX.CorrectionStds;

if strcmp(Plot,'yes')
    fold = "fig";

    if ~exist(fold, 'dir')
        mkdir(fold)
    end
    subfold = sprintf('./fig/%s',Mode);
    if ~exist(subfold, 'dir')
        mkdir(subfold)
    end
    set(0, 'DefaultFigureVisible', 'off');

    clf
    f = figure('Visible','off','Position',[1 50 1200 900]);
    tiledlayout("flow","TileSpacing","compact"); mx = max(dfs(:,1));
    for i = 1:length(udfs3)
        nexttile
        dy = DFC(i).X(:,1); dx = DFC(i).X(:,2);
        dxy = sortrows([dx dy],2);
        plot(dxy(:,1),dxy(:,2),'-+k');
        ylim([0 mx]);
        title(sprintf('C%s',string(udfs3(i))))
        set(gca,'YDir','reverse')
    end
    saveas(f,sprintf('./%s/%s/AllComponents_%s.png',fold,Mode,Mode))
    set(0, 'DefaultFigureVisible', 'on');

    if strcmp(Plot,'no'); end

[AQ.Function] = deal('assignQ');

end

