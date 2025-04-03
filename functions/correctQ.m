function [QX] = correctQ(AX,IX,varargin)
% correctQ - Applies outlier analysis and corrections to samples.
%
% This function performs compound-specific or bulk drift corrections, applies 
% size corrections, and normalizes scale using either a global regression or 
% two-point method. Additionally, outlier removal can be applied before 
% corrections, with multiple detection methods available, as set in infoQ.
%
% Inputs:
%   AX              - Output from findQ [struct]
%   IX              - Output from infoQ [struct]
%                     Defaults will be used for StandardID, PolyOrder, Mode, and OutlierMethod if not specified.
%   Plot            - Option to generate output figures into a directory ('yes' or 'no')
%
% Outputs:
%   QX              - Struct with fields:
%                     QX.Sample         - Corrected sample data
%                     QX.Standard       - Corrected standard data
%                     QX.Outlier        - Table of identified outliers (if applicable)
%                     QX.CorrectionStds - Table of drift, size, and scale correction standard deviations
%                     QX.FileName       - Processed file name
%
% Example Usage:
%   QX = correctQ(AX, IX, 'Plot', 'yes');
%
%   % Runs correction using the outlier method specified in infoQ.
%
%   QX = correctQ(AX, IX, 'Plot', 'no');
%
%   % Runs correction without generating plots.
%
% Notes:
%   - Outlier detection method is now set in **infoQ**, not correctQ.
%   - Available outlier detection methods (set in infoQ):
%       - 'mad'    - Median Absolute Deviation (Default, robust)
%       - 'zscore' - Z-score (removes values >3 standard deviations)
%       - 'iqr'    - Interquartile Range (removes values outside 1.5×IQR range)
%       - 'grubbs' - Grubbs' test (detects single extreme outliers in small datasets)
%       - 'none'   - No outlier removal
%   - Outlier detection is performed **before** applying drift, size, and scale corrections.
%   - The default outlier detection method is 'mad' (Median Absolute Deviation).
%   - To disable outlier removal, set 'OutlierMethod', 'none' in infoQ.

defPlot = 'yes';

expPlot = {'yes','no'};

p = inputParser;
validAX = @(x) isstruct(x);
validIX = @(x) isstruct(x);
validPlot = @(x) any(validatestring(x,expPlot));

addRequired(p,'AX',validAX)
addRequired(p,'IX',validIX)
addParameter(p,'Plot',defPlot,validPlot)

parse(p,AX,IX,varargin{:})

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

AX = p.Results.AX;
IX = p.Results.IX;
Plot = char(p.Results.Plot);

GX = AX;
std_id = IX.StandardID;
pln = IX.PolyOrder;
Mode = IX.Mode;
dcm = IX.DriftCorrection;
scm = IX.ScaleCorrection;
scl2 = IX.Scale2PTComp;
otm = IX.OutlierMethod;
otmtr = IX.OutlierThreshold;

QX = struct();

for j = 1:length(GX)
    G = GX(j).Analysis;
    lab = G.Identifier;

    temp = char(lab);
    idi = string(temp(:,1:2)) == std_id;

    Gstd = G(idi,:);
    Gsmp = G(~idi,:);

    rstd = Gstd.KnownValue;
    dstd = Gstd.RawDeltaValue;
    astd = Gstd.Analysis;
    cstd = Gstd.Component;
    zstd = Gstd.PeakAmplitude;
    pstd = Gstd.PeakArea;
    diffstd = dstd - rstd;

    rsmp = Gsmp.KnownValue;
    dsmp = Gsmp.RawDeltaValue;
    asmp = Gsmp.Analysis;
    csmp = Gsmp.Component;
    zsmp = Gsmp.PeakAmplitude;
    psmp = Gsmp.PeakArea;

    %%% Remove Outliers for Sample Data (if enabled)
    if strcmp(otm, 'mad')
        % Compute MAD-based threshold for sample data
        med_val = median(dsmp, 'omitnan');
        mad_val = median(abs(dsmp - med_val), 'omitnan');
        madthr = otmtr;
        dsmp_out = dsmp;

        % Identify outliers
        outliers = abs(dsmp - med_val) > madthr * mad_val;

        % Remove outliers (set to NaN)
        dsmp(outliers) = NaN;

        % Store outlier indices and values
        outindx = find(outliers);
        outval = dsmp_out(outliers);
        outtab = table(outindx, outval, 'VariableNames', {'Index', 'Value'});
    elseif strcmp(otm, 'zscore')
        z_scores = (dsmp - mean(dsmp, 'omitnan')) ./ std(dsmp, 'omitnan');
        z_thresh = otmtr; % 3 standard deviations
        outliers = abs(z_scores) > z_thresh;
        dsmp_out = dsmp;
        dsmp(outliers) = NaN;
        outindx = find(outliers);
        outval = dsmp_out(outliers);
        outtab = table(outindx, outval, 'VariableNames', {'Index', 'Value'});
    elseif strcmp(otm, 'iqr')
        Q1 = prctile(dsmp, 25); % First quartile (25th percentile)
        Q3 = prctile(dsmp, 75); % Third quartile (75th percentile)
        IQR_val = Q3 - Q1;
        lower_bound = Q1 - 1.5 * IQR_val;
        upper_bound = Q3 + 1.5 * IQR_val;

        outliers = (dsmp < lower_bound) | (dsmp > upper_bound);
        dsmp_out = dsmp;
        dsmp(outliers) = NaN;
        outindx = find(outliers);
        outval = dsmp_out(outliers);
        outtab = table(outindx, outval, 'VariableNames', {'Index', 'Value'});
    elseif strcmp(otm, 'grubbs')
        alpha = 0.05; % Significance level
        n = length(dsmp);
        t_crit = tinv(1 - alpha / (2 * n), n - 2); % Critical t-value
        G_thresh = ((n - 1) / sqrt(n)) * sqrt(t_crit^2 / (n - 2 + t_crit^2));

        mean_dsmp = mean(dsmp, 'omitnan');
        std_dsmp = std(dsmp, 'omitnan');
        G_scores = abs(dsmp - mean_dsmp) ./ std_dsmp;
        dsmp_out = dsmp;
        outliers = G_scores > G_thresh;

        dsmp(outliers) = NaN;
        outindx = find(outliers);
        outval = dsmp_out(outliers);
        outtab = table(outindx, outval, 'VariableNames', {'Index', 'Value'});
    elseif strcmp(otm, 'none')
        % No outliers removed
        outtab = table([], [], 'VariableNames', {'Index', 'Value'});
    end

    cuniq = unique(cstd); cuniqs = unique(csmp);
    if strcmp(dcm,'CompoundSpecific')
        % drift correction (Component-specific)
        polyn = pln;
        dstd_drift = zeros([length(dstd) 1]);
        dsmp_drift = zeros([length(dsmp) 1]);
        p_drift = zeros([polyn+1 length(cuniq)]);
        for i = 1:length(cuniq)
            ci = cstd == cuniq(i);
            di = dstd(ci);
            ai = astd(ci);
            dif = diffstd(ci);
            [pdi,~] = polyfit(ai,dif,polyn);
            dci = polyval(pdi,ai);
            ddi = di - dci;
            dstd_drift(ci) = ddi;
            p_drift(:,i) = pdi';

            ck = csmp == cuniq(i);
            dk = dsmp(ck);
            ak = asmp(ck);
            dck = polyval(pdi,ak);
            ddk = dk - dck;
            dsmp_drift(ck) = ddk;
        end
        cq = cuniqs > max(cuniq) | cuniqs < min(cuniq);
        cqn = cuniqs(cq); %%%%%%%%%%%% *apply C30 correction for C31-35
        for i = 1:length(cqn)
            if cqn(i) > max(cuniq)
                cm = csmp == cqn(i);
                dm = dsmp(cm);
                am = asmp(cm);
                pdm = p_drift(:,end);
                dcmr = polyval(pdm,am);
                ddm = dm - dcmr;
                dsmp_drift(cm) = ddm;
            elseif cqn(i) < min(cuniq)
                cm = csmp == cqn(i);
                dm = dsmp(cm);
                am = asmp(cm);
                pdm = p_drift(:,1);
                dcmr = polyval(pdm,am);
                ddm = dm - dcmr;
                dsmp_drift(cm) = ddm;
            end
        end
    end
    if strcmp(dcm,'Global')
        [pdi,~] = polyfit(astd,diffstd,pln);
        dci = polyval(pdi,astd);
        ddi = dstd - dci;
        dstd_drift = ddi;

        dck = polyval(pdi,asmp);
        ddk = dsmp - dck;
        dsmp_drift = ddk;
    end
    dsmp_drift(dsmp_drift == 0) = NaN;
    res_corr_smp1 = dsmp_drift - dsmp;
    res_corr_dstd = dstd_drift - dstd;
    res_iub1 = dstd_drift - rstd;
    resnan = ~isnan(res_iub1);
    drift_std = std(res_corr_dstd(~isnan(res_corr_dstd)));

    % size correction
    [ps,~] = polyfit(zstd(resnan),res_iub1(resnan),1);
    zidf = polyval(ps,zsmp);
    dsmp_size = dsmp_drift - zidf;
    res_corr_size = dsmp_size - dsmp_drift;
    size_std = std(res_corr_size(~isnan(res_corr_size)));

    % scale correction (expansion)
    if strcmp(scm,'Regression')
        valid_std = ~any(isnan([dstd_drift, rstd, pstd]), 2); % Ensure valid data for standards
        dstd_valid = dstd_drift(valid_std);
        rstd_valid = rstd(valid_std);
        pstd_valid = pstd(valid_std);
        wnorm = pstd_valid / max(pstd_valid); % Normalized weights
        wdiag = diag(wnorm); % Diagonal weight matrix
        dr = [dstd_valid(:), ones(length(dstd_valid), 1)]; % Matrix of measured standards
        fitpw = (dr' * wdiag * dr) \ (dr' * wdiag * rstd_valid);
        rr = fitpw(1); % Stretch factor (slope)
        offw = fitpw(2); % Offset (intercept)
        dsmp_scale = (dsmp_size - offw) ./ rr;
        res_corr_scl = dsmp_scale - dsmp_size;
    end
    if strcmp(scm,'TwoPoint')
        nct1 = scl2(1); nct2 = scl2(2);
        dstd_range = abs(mean(dstd(cstd == nct2)) - mean(dstd(cstd == nct1)));
        rep_range = abs(mean(rstd(cstd == nct2)) - mean(rstd(cstd == nct1)));
        rr = rep_range/dstd_range;  %%% stretching factor
        dsmp_r = dsmp_size*rr;
        dstd_r = dstd_drift*rr;
        ds1 = mean(mean(rstd(cstd == nct1) - dstd_r(cstd == nct1)));
        ds2 = mean(mean(rstd(cstd == nct2) - dstd_r(cstd == nct2)));
        ds = mean([ds1 ds2]);
        dsmp_scale = dsmp_r - ds;
        res_corr_scl = dsmp_scale - dsmp_size;
    end
    scale_std = std(res_corr_scl(~isnan(res_corr_scl)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Qsmp = Gsmp; Qstd = Gstd;
    Qsmp.CorrectedDeltaValue = dsmp_scale;
    Qstd.CorrectedDeltaValue = dstd_drift;

    cb = [57 106 150]./255; ms = 50; bcol = 'parula';

    if strcmp(Plot,'yes')
        set(0, 'DefaultFigureVisible', 'off');
        fold = "./fig";
        if ~exist(fold, 'dir')
            mkdir(fold)
        end
        subfold = sprintf('./fig/%s',Mode);
        if ~exist(subfold, 'dir')
            mkdir(subfold)
        end

        clf
        f = figure('Visible','off','Position',[210 55 970 777]);
        subplot(3,1,1)
        scatter(asmp,res_corr_smp1,ms/2,csmp,'filled','o');
        ylabel('Drift Adjustment (‰)');
        title(sprintf('File: "%s"',GX(j).FileName),'FontWeight','normal','Interpreter','none')
        cbr = colorbar; colormap(bcol);  title(cbr,'nC #');
        subplot(3,1,2)
        scatter(asmp,res_corr_size,ms/2,csmp,'filled','o');
        ylabel('Size Adjustment (‰)');
        cbr = colorbar; colormap(bcol);
        subplot(3,1,3)
        scatter(asmp,res_corr_scl,ms/2,csmp,'filled','o');
        ylabel('Scale Adjustment (‰)');
        cbr = colorbar; colormap(bcol);
        xlabel('Analysis');
        saveas(f,sprintf('./%s/%s/Corr%d_%s.png',fold,Mode,j,Mode))
        if strcmp(otm,'mad')
            f = figure('Visible','off','Position',[210 55 970 777]);
            subplot(2, 1, 1);
            plot(dsmp_out, '-xk', 'LineWidth', 0.9); hold on
            plot(outindx, outval, 'ro', 'MarkerSize', 8,'LineWidth', 1.6); hold off % Plot outliers
            title('Original Data with Outliers');
            text(0.01,0.05,sprintf('File: "%s"',GX(j).FileName),'FontWeight', ...
                'normal','Interpreter','none','Units','normalized')
            xlabel('Index');
            ylabel('Value');
            subplot(2, 1, 2);
            plot(dsmp_scale, '-xk', 'LineWidth', 0.9);
            title('Data After Outlier Removal');
            xlabel('Index');
            ylabel('Value');
            saveas(f,sprintf('./%s/%s/Outliers%d_%s.png',fold,Mode,j,Mode))
        end

        set(0, 'DefaultFigureVisible', 'on');
    end
    if strcmp(Plot,'no'); end

    QX(j).Sample = Qsmp;
    QX(j).Standard = Qstd;
    QX(j).Outlier = outtab;
    QX(j).CorrectionStds = table(drift_std,size_std,scale_std);
    QX(j).FileName = GX(j).FileName;

end

[QX.Mode] = deal(Mode);
[QX.Function] = deal('correctQ');

if strcmp(Plot,'yes')
    set(0, 'DefaultFigureVisible', 'off');

    f = figure('Visible','off','Position',[210 55 970 777]);clf
    avgpa = zeros(size(cuniqs));
    for i = 1:numel(cuniqs)
        cl = cuniqs(i);
        avgpa(i) = mean(Qsmp.PeakArea(Qsmp.Component == cl));
    end
    bar(cuniqs,avgpa,'FaceColor','#022688');
    xlabel('Compound');
    ylabel('Average Peak Area (Vs)');
    grid minor;
    saveas(f,sprintf('./%s/%s/AveragePeakArea_%s.png',fold,Mode,Mode))


    f = figure('Visible','off','Position',[210 55 970 777]);clf
    plot(Qstd.RawDeltaValue,Qstd.KnownValue,'.','Color','#022688','MarkerSize',12);
    xlabel('\delta (Measured)');
    ylabel('\delta (Known)');
    text(0.05,0.95,sprintf('Outliers Removed: %d',length(outval)),'Units','normalized');
    saveas(f,sprintf('./%s/%s/Offset_%s.png',fold,Mode,Mode))

    set(0, 'DefaultFigureVisible', 'on');

end
if strcmp(Plot,'no'); end


end

