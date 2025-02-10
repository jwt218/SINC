function [WT] = weightQ(AQ,IX,varargin)
% weightQ - Determines weighted average of selected compounds using peak
%           areas.
% Inputs:
%   AQ              - Output from assignQ [struct]
%   IX              - Output from infoQ [struct]
%                     Defaults will be used for Mode, Weight, OutputFile, Comp,
%                     and AddComp if not specified.
%   Error             - Output from errorQ [struct]
%   Plot            - Option to generate output figures into a directory ('yes' or 'no')
%   RemoveFlags     - Option to remove values above uncertainty threshold ('yes' or 'no')
% Outputs: [WT] is a struct with fields Data.
%
% Example Usage:
%   WT = weightQ(AQ, IX, 'Error', UCX, 'Plot', 'yes', 'RemoveFlags', 'yes');

defUCX = struct('TotalUncertainty',1);
defPlot = 'yes';
defRemoveFlags = 'no';

expPlot = {'yes','no'};
expRemoveFlags = {'yes','no'};

p = inputParser;
validAQ = @(x) isstruct(x);
validUCX = @(x) isstruct(x);
validIX = @(x) isstruct(x);
validPlot = @(x) any(validatestring(x,expPlot));
validRemoveFlags = @(x) any(validatestring(x,expRemoveFlags));

addRequired(p,'AQ',validAQ)
addRequired(p,'IX',validIX)
addParameter(p,'Error',defUCX,validUCX)
addParameter(p,'Plot',defPlot,validPlot)
addParameter(p,'RemoveFlags',defRemoveFlags,validRemoveFlags)

parse(p,AQ,IX,varargin{:})

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

AQ = p.Results.AQ;
IX = p.Results.IX;
UCX = p.Results.Error;
Plot = char(p.Results.Plot);
RemoveFlags = p.Results.RemoveFlags;

Mode = IX.Mode;
cm = IX.WeightedComp;
cmp = [IX.Compounds;IX.AdditionalComp];
acmp = IX.AdditionalComp;
age = AQ.Age;

exceedingValues = cm(cm > max(cmp)); % Values in A greater than max(B)
belowValues = cm(cm < min(cmp));     % Values in A lower than min(B)

% Check for violations
if ~isempty(exceedingValues) || ~isempty(belowValues)
    % Construct error message
    errMsg = "Error: Some weight compounds are outside the range of Comp + AddComp. ";
    
    if ~isempty(exceedingValues)
        errMsg = errMsg + sprintf("Values exceeding range (C%d): C%s.", max(cmp), num2str(exceedingValues));
    end
    
    if ~isempty(belowValues)
        errMsg = errMsg + sprintf("Values below range (C%d): C%s.", min(cmp), num2str(belowValues));
    end
    
    error(errMsg); % Throw error with formatted message
end


if strcmp(RemoveFlags,'no')
    K = AQ.MeanDelta;
    E = AQ.SDDelta;
    if isnumeric(UCX.TotalUncertainty)
        I = array2table(zeros(size(K)));
    else
        I = UCX.TotalUncertainty;
    end
    S = AQ.PeakArea;
    ctab = table2array(K(:,3:end));
elseif strcmp(RemoveFlags,'yes')
    K = UCX.FlaggedRemovedDelta;
    E = AQ.SDDelta;
    I = UCX.FlaggedRemovedUncertainty;
    S = AQ.PeakArea;
    ctab = table2array(K(:,2:end));
end

ctun = K.ID;

act = age;
stab = table2array(I(:,2:end));
ptab = table2array(S(:,3:end));
etab = table2array(E(:,3:end));
wk = ctab(:,ismember(cmp,cm));
sk = stab(:,ismember(cmp,cm));
pk = ptab(:,ismember(cmp,cm));
ek = etab(:,ismember(cmp,cm));
wknan = ~isnan(wk); sknan = ~isnan(sk);


if any(ismember(cm,acmp))
    eki = ismember(cm,acmp);
    pek = ek(:,eki);
    sk(:,eki) = pek.*2;
end

wm = zeros([length(ctun) 1]); ws = zeros([length(ctun) 1]);
for i = 1:length(ctun)
    wmj = sum(wk(i,wknan(i,:)).*pk(i,wknan(i,:)))/sum(pk(i,wknan(i,:)));
    wsj = sum(sk(i,sknan(i,:)).*pk(i,sknan(i,:)))/sum(pk(i,sknan(i,:)));
    wm(i) = wmj; ws(i) = wsj;
end

wtvn = {'ID','Age','Weighted mean','Weighted SD',...
    sprintf('Delta %s',char(strjoin(string(cm)))),...
    sprintf('Uncertainty %s',char(strjoin(string(cm)))),...
    sprintf('Weights %s',char(strjoin(string(cm))))};
WTu = table(ctun,act,wm,ws,wk,sk,pk,'VariableNames',wtvn);
WTx = sortrows(WTu,"Age");

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

    for i = 1:length(cm)
        f = figure('Visible','off','Position',[135 67 754 754]);
        cc = '#154360'; cc2 = '#922b21';
        xii = isfinite(wk(:,i)); ddx = wk(:,i); iix = sk(:,i);
        y = WTx.Age(xii); x = ddx(xii); wy = iix(xii);
        if any(cm(i) == acmp)
            plot(x,y,'-o','Color',cc2,'LineWidth',1.3,'MarkerFaceColor','w'); hold on
            errorbar(x,y,wy/2,'horizontal','LineStyle','none','Color',cc2,'Marker','.')
            hold off
        else
            plot(x,y,'-o','Color',cc,'LineWidth',1.3,'MarkerFaceColor','w'); hold on
            errorbar(x,y,wy/2,'horizontal','LineStyle','none','Color',cc,'Marker','.')
            hold off
        end
        xlabel('\delta (‰)'); ylabel('Age'); box on;
        legend('\delta','\sigma','box','off','Location','best')
        text(0.01,0.05,sprintf('nC#: %s',string(cm(i))),'Units','normalized')
        set(gca,'YDir','reverse','Position',[.08 .06 .9 .91]);
        fontsize(f, 12, "points")
        saveas(f,sprintf('./%s/%s/DeltaC%d_%s.png',fold,Mode,cm(i),Mode))
    end

    clf
    f = figure('Visible','off','Position',[135 67 754 754]);
    cc = '#17202a';
    xii = isfinite(WTx.("Weighted mean"));
    y = WTx.Age(xii); x = WTx.("Weighted mean")(xii); wy = WTx.("Weighted SD")(xii);
    plot(x,y,'-o','Color',cc,'LineWidth',1.3,'MarkerFaceColor','w'); hold on
    errorbar(x,y,wy/2,'horizontal','LineStyle','none','Color',cc,'Marker','.')
    hold off
    xlabel('Weighted \delta (‰)'); ylabel('Age'); box on;
    legend('\delta','\sigma','box','off','Location','best')
    text(0.01,0.05,sprintf('Weight Components: %s',strjoin(string(cm))),'Units','normalized')
    set(gca,'YDir','reverse','Position',[.08 .06 .9 .91]);
    fontsize(f, 12, "points")
    saveas(f,sprintf('./%s/%s/WeightedMean_%s.png',fold,Mode,Mode))
end
set(0, 'DefaultFigureVisible', 'off');

if strcmp(Plot,'no'); end

WT = struct();
WT.Data = WTx;
[WT.Function] = deal('weightQ');

end
