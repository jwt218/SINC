function [CX] = combineQ(QX,varargin)
% combineQ - Merges replicate analyses using common sample identifiers.
%
% Inputs:
%   QX       - Output from correctQ [struct]
%   Plot     - Option to generate output figures into a directory ('yes' or 'no')
%
% Outputs:
%   CX       - Struct with fields:
%              CX.Sample - Merged sample data
%              CX.Standard - Merged standard data
%              CX.CorrectionStds - Standard deviations of corrections
%
% Example Usage:
%   CX = combineQ(QX, 'Plot', 'yes');

defPlot = 'yes';

expPlot = {'yes','no'};

p = inputParser;
validQX = @(x) isstruct(x);
validPlot = @(x) any(validatestring(x,expPlot));

addRequired(p,'QX',validQX)
addParameter(p,'Plot',defPlot,validPlot)

parse(p,QX,varargin{:})

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

QX = p.Results.QX;
Plot = char(p.Results.Plot);
Mode = QX(1).Mode;

CX = struct();
Q = table(); Qs = table(); QCS = table();
for j = 1:length(QX)
    Q = [Q;QX(j).Sample];
    Qs = [Qs;QX(j).Standard];
    QCS = [QCS;QX(j).CorrectionStds];
end
sn = Q.Identifier;
sns = Qs.Identifier;

ll = length(Q.Analysis); lls = length(Qs.Analysis);
anys = zeros([ll 1]); anyss = zeros([lls 1]); qanys = ones([ll 1]); qanyss = ones([lls 1]);
strun = string(unique(sn(cellfun('isclass',sn,'char')),'stable'));
struns = string(unique(sns(cellfun('isclass',sns,'char')),'stable'));
for i = 1:length(strun)
    sni = strun(i);
    ki = contains(sn,sni);
    anys(ki) = i;
end
for i = 1:length(struns)
    sni = struns(i);
    ki = contains(sns,sni);
    anyss(ki) = i;
end
dqn = [0;diff(Q.Analysis)];
for i = 1:length(dqn)
    if i == 1
        qanys(i) = 1;
    elseif dqn(i) < 0
        qanys(i) = qanys(i-1) + 1;
    else
        qanys(i) = qanys(i-1);
    end
end
dqa = [0;diff(Qs.Analysis)];
for i = 1:length(dqa)
    if i == 1
        qanyss(i) = 1;
    elseif dqa(i) < 0
        qanyss(i) = qanyss(i-1) + 1;
    else
        qanyss(i) = qanyss(i-1);
    end
end

Q.SN = anys; Qs.SN = anyss;
Q.SEQ = qanys; Qs.SEQ = qanyss;
CX.Sample = Q;
CX.Standard = Qs;
CX.CorrectionStds = QCS;


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

    clf; f = figure('Visible','off','Position',[1 49 1920 955]);
    c = Q.component; ms = 100; cn = 29; comp = sort(unique(c));
    colormap(lines(3)); cf = find(c == cn);
    xv = xline(unique(anys),'--','Alpha',0.1); hold on
    s = scatter(anys,Q.d_final,ms,c,'filled','s');
    cb = colorbar; title(cb,'nC #'); cb.Ticks = comp(:);
    p = plot(anys(cf),Q.d_final(cf),'ok','MarkerSize',ms/8); hold off
    xlim([0 max(anys)+1])
    xlabel('Sample ID'); ylabel('\delta (‰; Corrected Values)')
    xticks(unique(anys,'stable'));
    xticklabels(unique(Q.Identifier1,'stable'));
    set(gca,'TickLabelInterpreter','none')
    legend(p,sprintf('nC%d, %d analyses',cn,length(strun)),'box','off','Location','northoutside')
    saveas(f,sprintf('./%s/%s/AllDelta_%s.png',fold,Mode,Mode))

    clf; f = figure('Visible','off','Position',[1 49 1920 955]);
    c = Q.component; ms = 60; cn = 29; comp = sort(unique(c));
    colormap(lines(3)); cf = find(c == cn);
    res = Q.("d_raw") - Q.d_final;
    xv = xline(unique(anys),'--','Alpha',0.1); hold on
    s = scatter(anys,res,ms,c,'filled','o');
    cb = colorbar; title(cb,'nC #'); cb.Ticks = comp(:); hold off
    xlim([0 max(anys)+1])
    xlabel('Sample ID'); ylabel('\delta (‰; Correction Residual)')
    xticks(unique(anys,'stable'));
    xticklabels(unique(Q.Identifier1,'stable'));
    set(gca,'TickLabelInterpreter','none')
    saveas(f,sprintf('./%s/%s/SampleResidual_%s.png',fold,Mode,Mode))

    clf; f = figure('Visible','off','Position',[1 49 1920 955]);
    c = Q.component; ms = 60; cn = 29; comp = sort(unique(c));
    colormap(lines(3)); cf = find(c == cn);
    res = Q.("d_raw") - Q.d_final;
    xv = xline(unique(qanys),'--','Alpha',0.1); hold on
    s = scatter(qanys,res,ms,c,'filled','o');
    cb = colorbar; title(cb,'nC #'); cb.Ticks = comp(:); hold off
    xlim([0 max(qanys)+1])
    xlabel('Sequence'); ylabel('\delta (‰; Correction Residual)')
    xticks(unique(qanys,'stable'));
    xticklabels(unique(Q.SEQ,'stable'));
    set(gca,'TickLabelInterpreter','none')
    saveas(f,sprintf('./%s/%s/SequenceResidual_%s.png',fold,Mode,Mode))

    clf; f = figure('Visible','off','Position',[30 150 1070 650]);
    c = Qs.component; ms = 60; cn = 29; comp = sort(unique(c));
    colormap(lines(3)); cf = find(c == cn);
    peaka = Qs.Amplitude; off = Qs.reported_IUB - Qs.d_raw;
    s = scatter(peaka,off,ms,c,'filled','s');
    cb = colorbar; title(cb,'nC #'); cb.Ticks = comp(:);
    xlabel('Amplitude'); ylabel('\delta (‰; Against Reference Material)')
    set(gca,'TickLabelInterpreter','none')
    saveas(f,sprintf('./%s/%s/StandardAmp_%s.png',fold,Mode,Mode))

    clf; f = figure('Visible','off','Position',[30 150 1070 650]);
    c = Q.component; ms = 60; cn = 29; comp = sort(unique(c));
    colormap(lines(3)); cf = find(c == cn);
    peaka = Q.Amplitude; off = Q.d_final - Q.reported_IUB;
    s = scatter(peaka,off,ms,c,'filled','s');
    cb = colorbar; title(cb,'nC #'); cb.Ticks = comp(:);
    xlabel('Amplitude'); ylabel('\delta (‰; Against Reference Material)')
    set(gca,'TickLabelInterpreter','none')
    saveas(f,sprintf('./%s/%s/SampleAmp_%s.png',fold,Mode,Mode))
    set(0, 'DefaultFigureVisible', 'on');

end
if strcmp(Plot,'no'); end
[CX.Function] = deal('combineQ');

end



