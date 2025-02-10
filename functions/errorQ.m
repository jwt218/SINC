function [UCX] = errorQ(AQ,AX,IX,varargin)
% errorQ - Calculates propagated uncertainty for each sample.
% Inputs:
%   AQ       - Output from assignQ [struct]
%   AX       - Output from findQ [struct]
%   IX       - Output from infoQ [struct]
%              Defaults will be used for SampleIDWin, Comp, AddComp,
%              ErrorThreshold, and Mode if not specified.
%   Plot     - Option to generate output figures into a directory ('yes' or 'no')
% Outputs: [UCX] is a struct with fields UCX.Delta, UCX.TotalUncertainty,
%          UCX.RefgasUncertainty, UCX.SampleUncertainty, UCX.FlaggedUncertainty,
%          UCX.FlaggedRemovedDelta, UCX.FlaggedRemovedUncertainty,
%          UCX.FlaggedSamples, and UCX.Covariance;
%
% Example Usage:
%   UCX = errorQ(AQ, AX, IX, 'Plot', 'yes');


defPlot = 'yes';

expPlot = {'yes','no'};

p = inputParser;
validCX = @(x) isstruct(x);
validAX = @(x) isstruct(x);
validIX = @(x) isstruct(x);
validPlot = @(x) any(validatestring(x,expPlot));

addRequired(p,'AQ',validCX)
addRequired(p,'AX',validAX)
addRequired(p,'IX',validIX)
addParameter(p,'Plot',defPlot,validPlot)

parse(p,AQ,AX,IX,varargin{:})

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

AQ = p.Results.AQ;
AX = p.Results.AX;
IX = p.Results.IX;
Plot = char(p.Results.Plot);
Mode = IX.Mode;

corstd = mean(table2array(AQ.CorrectionStds),1,'omitmissing');

tid = IX.SampleMatchLen;
A = table(); As = table();
for j = 1:length(AX)
    A = [A;AX(j).Analysis];
    As = [As;AX(j).Refgas];
end
L = AQ.Sample; Ls = AQ.Standard; Lgx = As; 
cmp = [IX.Compounds;IX.AdditionalComp];
thold = IX.UncertaintyThreshold;

dsmp = L.CorrectedDeltaValue;
csmp = L.Component;
rsmp = L.KnownValue;
repsd = L.KnownSD;
dref = Lgx.RawDeltaValue;
idref = Lgx.Identifier;

idsmp = char(L.Identifier);
idu = unique(string(idsmp(:,1:3)),'stable');
idx = string(idsmp(:,1:3));
idc = char(idx);
idrc = char(idref);
idrs = string(idrc(:,1:3));

age = AQ.Age;
alab = AQ.MeanDelta.ID;
lalb = char(idsmp);
ctar = lalb(:,1:tid);
stl = string(ctar);
[ctun,~] = unique(stl,'stable');
ai = ismember(alab,ctun);

dstd = Ls.CorrectedDeltaValue;
cstd = Ls.Component;

dxi = zeros([length(idu) length(cmp)]);
cxi = zeros([length(idu) length(cmp)]);
sxi = zeros([length(idu) length(cmp)]);
vxi = zeros([length(idu) length(cmp)]);
nxi = zeros([length(idu) length(cmp)]);
gxi = zeros([length(idu) length(cmp)]);
pxi = zeros([length(idu) length(cmp)]);
agi = zeros([length(idu) 1]);

for j = 1:length(idu)

    ki = string(idc) == string(idu(j));
    rxi = rsmp(ki);
    rxii = repsd(ki);
    drfx = dref(idrs == idu(j));
    dff = mean(drfx); if isnan(dff); dff = 0; end
    sff = std(drfx); if isnan(sff); sff = 0; end
    cmpki = csmp(ki); dki = dsmp(ki); Ca = [];

    for i = 1:numel(cmp)
        cd = find(cmpki == cmp(i));
        dkc = dki(cd);
        dm = mean(dkc);
        cmm = mean(cmpki(cd));
        smm = std(dkc,'omitnan');
        if cd < 2; smm = NaN; end

        rmm = mean(rxi(cd));
        rss = mean(rxii(cd));
        dss = dstd(cstd == cmm);
        ds = mean(dss);
        sss = sqrt(std(dss)^2 + corstd(1)^2 + corstd(2)^2 + corstd(3)^2);

        % Normalize the sample and reference gas values separately
        normalized_dm = dm - mean(dkc); % Normalize sample mean
        normalized_dff = dff - mean(dff); % Normalize refgas mean

        % Calculate ssM using the normalized values
        term1 = ((normalized_dm + 1)^2) * (sff^2); % Contribution from reference gas
        term2 = ((normalized_dff + 1)^2) * (smm^2); % Contribution from sample replicates
        ssM = term1 + term2; % Total sample uncertainty (squared)
        srM = (((1/(ds + 1))^2) * (rss^2)) + ((-(((rmm) + 1)/(((ds + 1))^2))^2) * (sss^2));
        covsr = 2 * term1 * (((1/(ds + 1))^2) * (rss^2)) * sff^2;

        sigx2 = ssM + srM;
        pcov = 2*(covsr/sigx2); % should be less than 1
        Ca = [Ca; cmm dm smm ds sss sqrt(abs(ssM)) sqrt(abs(srM)) sqrt(abs(sigx2)) pcov];
    end

    cx = Ca(:,1); cxi(j,:) = cx';
    dx = Ca(:,2); dxi(j,:) = dx';
    sx = Ca(:,3); sxi(j,:) = sx';
    sm = Ca(:,6); vxi(j,:) = sm';
    sl = Ca(:,7); nxi(j,:) = sl';
    su = Ca(:,8); gxi(j,:) = su';
    pc = Ca(:,9); pxi(j,:) = pc';
    agi(j) = age(idu(j) == alab);

end

cxi = sortrows([agi cxi],1); cxi = cxi(:,2:end);
dxi = sortrows([agi dxi],1); dxi = dxi(:,2:end);
sxi = sortrows([agi sxi],1); sxi = sxi(:,2:end);
vxi = sortrows([agi vxi],1); vxi = vxi(:,2:end);
nxi = sortrows([agi nxi],1); nxi = nxi(:,2:end);
gxi = sortrows([agi gxi],1); gxi = gxi(:,2:end);
pxi = sortrows([agi pxi],1); pxi = pxi(:,2:end);

txi = gxi;
txi(gxi < thold) = NaN;
txi2 = gxi;
txi2(gxi >= thold) = NaN;
kgx = gxi >= thold;
kgt = find(any(kgx > 0,2));

xxi = dxi;
xxi(kgx) = NaN;

vng = {'string','double','string','string'};
tgx = table('Size',[length(kgt) length(vng)],'VariableTypes',vng);
tgx.Properties.VariableNames = {'ID','#','nC Evens','nC Odds'};
for i  = 1:length(kgt)
    idg = idu(kgt(i));
    ng = length(find(kgx(kgt(i),:) > 0));
    cng = cxi(kgt(i),kgx(kgt(i),:) > 0);
    even = cng(rem(cng,2) == 0);
    odd = cng(rem(cng,2) == 1);
    tgx.ID(i) = idg;
    tgx.("#")(i) = ng;
    tgx.("nC Evens")(i) = strjoin(string(even));
    tgx.("nC Odds")(i) = strjoin(string(odd));
end


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
    f = figure('Visible','off','Position',[51 110 820 590]);
    ntxi = txi(:); nnt = length(ntxi(~isnan(ntxi))); nnx = length(ntxi(isnan(ntxi)));
    h = histogram(gxi(:),round(length(gxi(:))/50),'EdgeColor','none','FaceColor','r');
    xlabel('Total Analytical \sigma (‰)'); ylabel('#');
    p = patch([thold max(gxi(:))+1000 max(gxi(:))+1000 thold],[0 0 max(h.Values)+1000 max(h.Values)+1000],'k');
    p.FaceAlpha = 0.3; p.EdgeColor = "none";
    xline(thold,'--k','LineWidth',1.0); xlim([0 max(gxi(:))]); ylim([0 max(h.Values)])
    text(thold,max(h.BinCounts),sprintf('n = %d',nnt),'VerticalAlignment','top','HorizontalAlignment','right','Rotation',90)
    text(thold,max(h.BinCounts),sprintf('n = %d',nnx),'VerticalAlignment','bottom','HorizontalAlignment','right','Rotation',90)
    saveas(f,sprintf('./%s/%s/ErrHist.png',fold,Mode))

    clf
    f = figure('Visible','off','Position',[51 110 820 590]);
    ntxi = txi(:); nnt = length(ntxi(~isnan(ntxi))); nnx = length(ntxi(isnan(ntxi)));
    plot(cxi,gxi,'.','MarkerSize',18);
    ax = gca; ax.XTick = cmp;
    xlim([min(cmp) max(cmp)]); ylim([0 max(gxi(:))]);
    yline(thold,'--k','LineWidth',1.0);
    xlabel('nC#'); ylabel('Total Analytical \sigma (‰)');
    p = patch([min(cmp) max(cmp) max(cmp) min(cmp)],[thold thold max(gxi(:))+1000 max(gxi(:))+1000],'k');
    p.FaceAlpha = 0.3; p.EdgeColor = "none";
    text(min(cmp),thold,sprintf('n = %d',nnt),'VerticalAlignment','bottom')
    text(min(cmp),thold,sprintf('n = %d',nnx),'VerticalAlignment','top')
    saveas(f,sprintf('./%s/%s/TotalErr.png',fold,Mode))

    clf
    f = figure('Visible','off','Position',[51 110 820 590]);
    ntxi = txi(:); nnt = length(ntxi(~isnan(ntxi))); nnx = length(ntxi(isnan(ntxi)));
    scatter(vxi(:),nxi(:),22,cxi(:),'filled','o'); hold on
    plot([0 thold],[thold 0],'--k'); hold off
    xlim([min(vxi(:)) max(vxi(:))]); ylim([min(nxi(:)) max(nxi(:))]);
    xlabel('Sample \sigma (‰)'); ylabel('Refgas \sigma (‰)');
    p = patch([thold max(vxi(:)) max(vxi(:)) thold],[min(nxi(:)) min(nxi(:)) thold thold],'k');
    p.FaceAlpha = 0.3; p.EdgeColor = "none";
    text(thold,max(nxi(:)),sprintf('n = %d',nnt),'VerticalAlignment','top','HorizontalAlignment','right','Rotation',90)
    text(thold,max(nxi(:)),sprintf('n = %d',nnx),'VerticalAlignment','bottom','HorizontalAlignment','right','Rotation',90)
    cbr = colorbar; colormap('jet');  title(cbr,'nC #');
    saveas(f,sprintf('./%s/%s/RefgasSampleErr.png',fold,Mode))
    set(0, 'DefaultFigureVisible', 'on');

end
if strcmp(Plot,'no'); end

dz = [length(ai) length(cmp)];
vi = string(cmp);
jri = ["ID";vi];
tia = find(ai == 1);
gri = NaN(dz); vri = NaN(dz); nri = NaN(dz);

for i = 1:length(tia); gri(tia(i),:) = gxi(i,:); end
for i = 1:length(tia); vri(tia(i),:) = vxi(i,:); end
for i = 1:length(tia); nri(tia(i),:) = nxi(i,:); end
gci = [repmat(9999, dz(1), 1) gri];
vci = [repmat(9999, dz(1), 1) vri];
nci = [repmat(9999, dz(1), 1) nri];
U = array2table(gci); % combined uncertainty
U.Properties.VariableNames = jri;
U.ID = alab;
V = array2table(vci); % refgas uncertainty
V.Properties.VariableNames = jri;
V.ID = alab;
N = array2table(nci); % sample uncertainty
N.Properties.VariableNames = jri;
N.ID = alab;

tri = NaN(dz); xri = NaN(dz); tri2 = NaN(dz); dri = NaN(dz); pri = NaN(dz);
for i = 1:length(tia); tri(tia(i),:) = txi(i,:); end
for i = 1:length(tia); xri(tia(i),:) = xxi(i,:); end
for i = 1:length(tia); tri2(tia(i),:) = txi2(i,:); end
for i = 1:length(tia); dri(tia(i),:) = dxi(i,:); end
for i = 1:length(tia); pri(tia(i),:) = pxi(i,:); end
tci = [repmat(9999, dz(1), 1) tri];
F = array2table(tci);
F.Properties.VariableNames = jri;
F.ID = alab;
xci = [repmat(9999, dz(1), 1) xri];
X = array2table(xci);
X.Properties.VariableNames = jri;
X.ID = alab;
tci2 = [repmat(9999, dz(1), 1) tri2];
F2 = array2table(tci2);
F2.Properties.VariableNames = jri;
F2.ID = alab;
dzi = [repmat(9999, dz(1), 1) dri];
X2 = array2table(dzi);
X2.Properties.VariableNames = jri;
X2.ID = alab;
pzi = [repmat(9999, dz(1), 1) pri];
P = array2table(pzi);
P.Properties.VariableNames = jri;
P.ID = alab;

UCX = struct();
UCX.Delta = X2;
UCX.TotalUncertainty = U;
UCX.RefgasUncertainty = N;
UCX.SampleUncertainty = V;
UCX.FlaggedUncertainty = F;
UCX.FlaggedRemovedDelta = X;
UCX.FlaggedRemovedUncertainty = F2;
UCX.FlaggedSamples = tgx;
UCX.Covariance = P;

[UCX.Function] = deal('errorQ');

end
