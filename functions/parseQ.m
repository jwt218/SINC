function [EX] = parseQ(FileNames,varargin)
% parseQ - Reformats Qtegra .csv export file to improve readability and
%          apply proper format for other IsoNQ functions.
% Inputs:
%   FileNames       - List of file paths to be parsed [Fx1 string array]
%   Mode            - Specify parsing for d13C or dD data ('C' or 'H') [char]
% Outputs: [EX] is a struct with fields EX.Data and EX.FileName.
%
% Example Usage:
%   fnames = ["sample1.csv", "sample2.csv"];
%   EX = parseQ(fnames, 'Mode', 'C');

defMode = 'C';
expMode = {'C','H'};

p = inputParser;
validFileNames = @(x) isstring(x);
validMode = @(x) any(validatestring(x,expMode));

addRequired(p,'FileNames',validFileNames)
addParameter(p,'Mode',defMode,validMode)

parse(p,FileNames,varargin{:})

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

fnames = p.Results.FileNames(:);
Mode = p.Results.Mode(:);

warning('off')
EX = struct();

for j = 1:length(fnames)
    R = readtable(fnames(j),'ReadVariableNames',true,'VariableNamingRule','preserve');

    if strcmp(Mode,'C')
        sl = table2array(R(:,"Sample List - Vial"));   % vial pos
        sn = string(table2array(R(:,"Sample List - Label")));   % label
        rt = table2array(R(:,"RetentionTime 45.00 m/z - Value"));  % time 45 (s)
        pka45 = table2array(R(:,"PeakRawArea 45.00 m/z - Value")); % peak area
        pkh45 = table2array(R(:,"PeakAmplitude 45.00 m/z - Value")); % peak amplitude
        pkn = table2array(R(:,"Sample List - Run")); % peak number
        rco2 = table2array(R(:,"Ratio 45C.O.O/44C.O.O - Value")); % ratio
        bgd = table2array(R(:,"Background 45.00 m/z - Value")); % background

        exc = strcmp('d13 - Value',R.Properties.VariableNames);
        if any(exc); draw = table2array(R(:,"d13 - Value")); end
        exc = strcmp('d13C - Value',R.Properties.VariableNames);
        if any(exc); draw = table2array(R(:,"d13C - Value")); end

    elseif strcmp(Mode,'H')
        sl = table2array(R(:,"Sample List - Vial"));   % vial pos
        sn = string(table2array(R(:,"Sample List - Label")));   % label
        rt = table2array(R(:,"RetentionTime 3.00 m/z - Value"));  % time 45 (s)
        pka45 = table2array(R(:,"PeakRawArea 3.00 m/z - Value")); % peak area
        pkh45 = table2array(R(:,"PeakAmplitude 3.00 m/z - Value")); % peak amplitude
        pkn = table2array(R(:,"Sample List - Run")); % peak number
        rco2 = table2array(R(:,"Ratio 3H.H/2H.H - Value")); % ratio
        bgd = table2array(R(:,"Background 3.00 m/z - Value")); % background

        exc = strcmp('H2 d2 - Value',R.Properties.VariableNames);
        if any(exc); draw = table2array(R(:,"H2 d2 - Value")); end
        exc = strcmp('H2 d2H - Value',R.Properties.VariableNames);
        if any(exc); draw = table2array(R(:,"H2 d2H - Value")); end

    end
    strun = string(unique(sn(cellfun('isclass',sn,'char')),'stable'));
    ls = length(strun);

    ll = length(sl);
    anys = zeros([ll 1]);
    for i = 1:length(strun)
        sni = strun(i);
        ki = contains(sn,sni);
        anys(ki) = i;
    end

    vns = {'Identifier','Analysis','RetentionTime', ...
        'RawDeltaValue','PeakArea','PeakAmplitude','PeakNumber','Background'};
    vts = {'string','double','double','double','double',...
        'double','double','double'};
    E = table('Size',[ll,8],'VariableTypes',vts);

    E.Properties.VariableNames = vns;
    cnum = length(E.Properties.VariableNames);
    E(:,2:cnum) = array2table(repmat(9999, ll, cnum-1));
    for i = 1:ls
        sni = strun(i);
        ki = find(contains(sn,sni));
        pkni = pkn(ki);
        dpk = diff(pkni); dpkf = find(sign(dpk) == -1);

        l3 = length(ki)/3;
        li1 = ki(1:dpkf(1));
        li2 = ki(dpkf(1)+1:dpkf(2));
        li3 = ki(dpkf(2)+1:end);


        sli = anys(li1);  % analysis
        id1 = string(sn(li1)); % identifier 1
        rts = rt(li1);  % ret. time
        d13 = draw(li1); % raw d13C values
        pkr = pka45(li1); % peak raw area 45
        pkh = pkh45(li1); % peak amplitude 45
        pnn = pkn(li1); % peak number
        bck = bgd(li1); % background

        varcheck = {sli, id1, rts, d13, pkr, pkh, pnn, bck};
        nanflag = cellfun(@(x) isnumeric(x) && all(isnan(x)), varcheck);

        if any(nanflag)
            if nanflag(1), sli = anys(li2); end
            if nanflag(2), id1 = string(sn(li2)); end
            if nanflag(3), rts = rt(li2); end
            if nanflag(4), d13 = draw(li2); end
            if nanflag(5), pkr = pka45(li2); end
            if nanflag(6), pkh = pkh45(li2); end
            if nanflag(7), pnn = pkn(li2); end
            if nanflag(8), bck = bgd(li2); end
        end

        varcheck = {sli, id1, rts, d13, pkr, pkh, pnn, bck};
        nanflag = cellfun(@(x) isnumeric(x) && all(isnan(x)), varcheck);

        if any(nanflag)
            if nanflag(1), sli = anys(li3); end
            if nanflag(2), id1 = string(sn(li3)); end
            if nanflag(3), rts = rt(li3); end
            if nanflag(4), d13 = draw(li3); end
            if nanflag(5), pkr = pka45(li3); end
            if nanflag(6), pkh = pkh45(li3); end
            if nanflag(7), pnn = pkn(li3); end
            if nanflag(8), bck = bgd(li3); end
        end

        En = table(id1,sli,rts,d13,pkr,pkh,pnn,bck);
        En.Properties.VariableNames = vns;
        E(ki(1:size(En,1)),2:cnum) = array2table([sli rts d13 pkr pkh pnn bck]);
        E{ki(1:size(En,1)),1} = id1;
    end

    E(E.(3) == 9999,:) = [];
    EX(j).Data = E;
    EX(j).FileName = fnames(j);

end

[EX.Function] = deal('parseQ');

end
