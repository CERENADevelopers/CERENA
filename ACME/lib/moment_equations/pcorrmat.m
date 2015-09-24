% function [pCorrMat,pCorrMatAll] = pcorrmat(System,M,covMat,corrMatAll,options)
function [pCorrMat,pCorrMatAll] = pcorrmat(varargin)
if nargin>3
    System = varargin{1};
    M = varargin{2};
    covMat = varargin{3};
    corrMatAll = varargin{4};
else
    error('At least four input arguments are required!')
end

options.visualization = 'on';
options.movieName = 'pcorrMovie1';
options.frameRate = 15;
if nargin>4
    options = setdefault(varargin{5},options);
end

ndata = 100000;
simData = zeros(ndata,System.state.number);
s = zeros(System.state.number,System.state.number,size(M,1));
invCorrMatAll = zeros(System.state.number,System.state.number,size(M,1));
pCorrMatAll = zeros(System.state.number,System.state.number,size(M,1));
mu = M(:,1:System.state.number);
for j = 2:size(M,1)
        invCorrMatAll(:,:,j) = inv(corrMatAll(:,:,j));
        pCorrMatAll(:,:,j) = - invCorrMatAll(:,:,j)./sqrt(diag(invCorrMatAll(:,:,j))*diag(invCorrMatAll(:,:,j))');
    if (~isreal(pCorrMatAll(:,:,j)))
        simData(:,:,j) = mvnrnd(mu(j,:),covMat(:,:,j),ndata);
        [s(:,:,j), lamcor, lamvar] = covshrinkKPM(simData(:,:,j));
        corrMatAll(:,:,j) = corrcov(s(:,:,j));
        invCorrMatAll(:,:,j) = inv(corrMatAll(:,:,j));
        pCorrMatAll(:,:,j) = - invCorrMatAll(:,:,j)./sqrt(diag(invCorrMatAll(:,:,j))*diag(invCorrMatAll(:,:,j))');
    end 
end

%%
pCorrMat = zeros(size(pCorrMatAll,1),size(pCorrMatAll,1));
[~,ind_max_pcorr] = max(abs(pCorrMatAll),[],3);
for i = 1:size(pCorrMatAll,1)
    for j = 1:size(pCorrMatAll,1)
        pCorrMat(i,j) = pCorrMatAll(i,j,ind_max_pcorr(i,j));
    end
end
fs = 16;
    figure('Position',[1 1 900 1000]);
    set(gcf,'PaperPositionMode','Auto')
    imagesc(pCorrMat);
    set(gca, 'XTick', 1:System.state.number); % center x-axis ticks on bins
    set(gca, 'YTick', 1:System.state.number); % center y-axis ticks on bins
    set(gca, 'XTickLabel', System.state.name, 'FontSize', fs); % set x-axis labels
    set(gca, 'YTickLabel', System.state.name, 'FontSize', fs); % set y-axis labels
    xticklabel_rotate(1:System.state.number,45,System.state.name,'interpreter','none')
    set(gca, 'FontSize', fs);
    title('Partial correlation map', 'FontSize', fs); % set title
    colormap('hot');
    colorbar;
    set(gca, 'FontSize', fs);
if strcmp(options.visualization,'on')
    
    %% Visualization - partial correlation
    hfig = figure('Position',[1 1 900 1000]);
    set(gcf,'PaperPositionMode','Auto')
    for i = 1:size(pCorrMatAll,3)
        imagesc(pCorrMatAll(:,:,i));
        set(gca, 'XTick', 1:System.state.number); % center x-axis ticks on bins
        set(gca, 'YTick', 1:System.state.number); % center y-axis ticks on bins
        set(gca, 'XTickLabel', System.state.name, 'FontSize', fs); % set x-axis labels
        set(gca, 'YTickLabel', System.state.name, 'FontSize', fs); % set y-axis labels
        set(gca, 'FontSize', fs);
        title('Partial correlation map', 'FontSize', fs); % set title
        colormap('jet');
        colorbar;
        set(gca, 'FontSize', fs);
        frame(i) = getframe(hfig);
    end
    
    %%
    if isfield(options,'movieName')
        vidobj = VideoWriter([options.movieName,'.avi']);
    else
        vidobj = VideoWriter('pCorrMovie1.avi');
    end
    vidobj.FrameRate = options.frameRate;
    open(vidobj)
    for i =1:length(frame)
        writeVideo(vidobj,frame(i));
    end
    close(vidobj);
end