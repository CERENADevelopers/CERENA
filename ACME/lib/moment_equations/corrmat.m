% function [corrMat,corrMatAll,covMat] = corrmat(System,M,options)
function [corrMat,corrMatAll,covMat] = corrmat(varargin)
if nargin>1
    System = varargin{1};
    M = varargin{2};
else
    error('At least two input arguments are required!')
end

options.visualization = 'on';
options.movieName = 'corrMovie1';
options.frameRate = 15;
if nargin>2
    options = setdefault(varargin{3},options);
end

ndata = 100000;
simData = zeros(ndata,System.state.number);
s = zeros(System.state.number,System.state.number,size(M,1));
mu = M(:,1:System.state.number);


switch System.name
    case 'MM'
        state_order = System.MM.sym.state.order;
    case 'MCM'
        state_order = System.CMM.state_moments.order;
    case {'RRE','LNA','EMRE','IOS','FSP'}
        state_order = System.state.order;
end
ind_covar = find(sum(state_order~=0,2)==2);
order_cov = state_order(ind_covar,end-1:end);
ind_var = find((state_order(:,end-1)==state_order(:,end)).*(sum(state_order~=0,2)==2));

M_cov = M(:,ind_covar);

corr_coeff = M_cov./sqrt(M(:,ind_var(order_cov(:,1))).*M(:,ind_var(order_cov(:,2))));
[~,ind_max_corr_coeff] = max(abs(corr_coeff));
for i = 1:size(corr_coeff,2)
    corr_coeff_max(i) = corr_coeff(ind_max_corr_coeff(i),i);
end
% creating the covariance matrix and mean vector from the time points with the maximum correlation
covMat = zeros(System.state.number,System.state.number,size(M,1));
corrMatAll = zeros(System.state.number,System.state.number,size(M,1));
% mu1 = zeros(1,System.state.number);
for j = 2:size(M,1)
    k = 1;
    for i = 1:System.state.number
        covMat(i,i:end,j) = M_cov(j,k:k+System.state.number-i);
        corrMatAll(i,i:end,j) = corr_coeff(j,k:k+System.state.number-i);
        k = k+System.state.number-i+1;
    end
    covMat(:,:,j) = squeeze(covMat(:,:,j)) + tril(squeeze(covMat(:,:,j))', -1);
    corrMatAll(:,:,j) = squeeze(corrMatAll(:,:,j)) + tril(squeeze(corrMatAll(:,:,j))', -1);
    if (~isreal(corrMatAll(:,:,j)))
        simData(:,:,j) = mvnrnd(mu(j,:),covMat(:,:,j),ndata);
        [s(:,:,j), lamcor, lamvar] = covshrinkKPM(simData(:,:,j));
        corrMatAll(:,:,j) = corrcov(s(:,:,j));
    end
end
% creating the correlation matrix
corrMat = zeros(System.state.number,System.state.number);
k = 1;
for i = 1:System.state.number
    corrMat(i,i:end) = corr_coeff_max(k:k+System.state.number-i);
    k = k+System.state.number-i+1;
end
%
corrMat = corrMat + tril(corrMat', -1);
%% Visualization - max correlation map
    fs = 16;
    figure('Position',[1 1 900 1000]);
    set(gcf,'PaperPositionMode','Auto')
    imagesc(corrMat);
    set(gca, 'XTick', 1:System.state.number); % center x-axis ticks on bins
    set(gca, 'YTick', 1:System.state.number); % center y-axis ticks on bins
    set(gca, 'XTickLabel', System.state.name, 'FontSize', fs); % set x-axis labels
    set(gca, 'YTickLabel', System.state.name, 'FontSize', fs); % set y-axis labels
    xticklabel_rotate(1:System.state.number,45,System.state.name,'interpreter','none')
    set(gca, 'FontSize', fs);
    title('Correlation map', 'FontSize', fs); % set title
    colormap('jet');
    colorbar;
    set(gca, 'FontSize', fs);
if strcmp(options.visualization,'on')
    
    %% Visualization - correlation through time
    hfig = figure('Position',[1 1 900 1000]);
    set(gcf,'PaperPositionMode','Auto')
    % set(gca,'NextPlot','replaceChildren');
    for i = 1:size(corrMatAll,3)
        imagesc(corrMatAll(:,:,i));
        set(gca, 'XTick', 1:System.state.number); % center x-axis ticks on bins
        set(gca, 'YTick', 1:System.state.number); % center y-axis ticks on bins
        set(gca, 'XTickLabel', System.state.name, 'FontSize', fs); % set x-axis labels
        set(gca, 'YTickLabel', System.state.name, 'FontSize', fs); % set y-axis labels
        set(gca, 'FontSize', fs);
        title('Correlation map', 'FontSize', fs); % set title
        colormap('jet');
        colorbar;
        set(gca, 'FontSize', fs);
        frame(i) = getframe(hfig);  
    end

    %% Movie
    if isfield(options,'movieName')
        vidobj = VideoWriter([options.movieName, '.avi']);
    else
        vidobj = VideoWriter('corrMovie1.avi');
    end
    vidobj.FrameRate = options.frameRate;
    open(vidobj)
    for i =1:length(frame)
        writeVideo(vidobj,frame(i));
    end
    close(vidobj);
    
end