function prv=cAMPsecretion_cellmigration_LJpotential

    close all
    Export_path = uigetdir;
    ImageName = 'Pattern';
    %%% for visualization %%%%
    scrsz = get(groot,'ScreenSize');
    figure('visible','off','position',[0 0 scrsz(3)/2 scrsz(4)],...
        'NumberTitle','off','Name','cAMP field, chemtaxis, LJ potential');

    
    %%% condition %%%
    cellnum = 200;
    width = 50;
    tspan = 10000;
    dt = 0.001;

    %%%%% variables and function sets %%%%%
    % 3 variables
    cAMPfield = zeros(width,width);

    R = 0;
    while min(min(R)) < 1
        cellx = randi([1,width],cellnum,1);
        celly = randi([1,width],cellnum,1);

        dx = repmat(cellx,[1 cellnum]);
        dx = dx' - dx;
        dx(abs(dx) > width/2) = -(width - dx(abs(dx) > width/2));
        dy = repmat(celly,[1 cellnum]);
        dy = dy' - dy;
        dy(abs(dy) > width/2) = -(width - dy(abs(dy) > width/2));
        R = sqrt(dx.^2 + dy.^2);
        R(1:(size(R,1)+1):size(R,1)^2) = NaN;
    end

    prv = cell(3,1);
    prv{1} = cAMPfield;
    prv{2} = cat(2,cellx,celly);
    prv{3} = prv{2}*0; % position at previous frame.

    eqns = cell(2,1);
    eqns{1} = @cAMP_secdegdif;
    eqns{2} = @cellmig;

    imagesc(prv{1});
    hold on
    plot(prv{2}(:,1),prv{2}(:,2),'ro')
    axis image;
    text(1,5,['T=0'],'FontSize',20,'Color','white');
    colorbar; %caxis([0,1]);
    pause;

    %%%%%%%% simulation %%%%%%%%%
    for time = 1:tspan;

        %if time > 500; prv{1}(50:55,50:55)=1; end

        nxt = RunKut4(eqns,prv,dt);
        nxt{2} = pdc_calc(nxt{2},width);
        nxt{3} = nxt{2} - prv{2}; % velocity calculation
        nxt{3}(abs(nxt{3}) > width/2) = -(width - nxt{3}(abs(nxt{3}) > width/2));

        if rem(time,100) == 0
            hold off
            imagesc(nxt{1});
            hold on
            plot(nxt{2}(:,1),nxt{2}(:,2),'ro')
            axis image;
            text(1,5,['T=',num2str(time)],'FontSize',20,'Color','white');
            colorbar; %caxis([0,1]);
            %pause(1/100);
            saveas(gcf,[Export_path,'/test',num2str(time,'%04d'),'.jpg']);
            %imwrite(nxt(:,:,1),[Export_path1,'/',ImageName,num2str(j,'%03d'),'.tif']);
        end
        prv = nxt;
    end

end

%%%%%%%% Runge Kutta 4 %%%%%%%%
function nxt = RunKut4(eqns,prv,dt)
    
    % eqns : cell array of the handles of differential equations.
    % prv  : cell array of the variables.
    % dt   : scalar value of dt.

    eqnsn = length(eqns);
    k1 = cellmul(prv,0);
    k2 = cellmul(prv,0);
    k3 = cellmul(prv,0);
    k4 = cellmul(prv,0);
    
    for i = 1:eqnsn
        if ~isempty(eqns{i})
            eqntemp = eqns{i};
            k1{i} = eqntemp(prv);
        end
    end

    for i = 1:eqnsn
        if ~isempty(eqns{i})
            eqntemp = eqns{i};
            k2{i} = eqntemp(celladd(prv,cellmul(k1,dt/2.0)));
        end
    end
    
    for i = 1:eqnsn
        if ~isempty(eqns{i})
            eqntemp = eqns{i};
            k3{i} = eqntemp(celladd(prv,cellmul(k2,dt/2.0)));
        end
    end
    
    for i = 1:eqnsn
        if ~isempty(eqns{i})
            eqntemp = eqns{i};
            k4{i} = eqntemp(celladd(prv,cellmul(k3,dt)));
        end
    end
    nxt = celladd(prv,cellmul(celladd(k1,k2,k3,k4),dt/6));
end

%%%%%%%% differential eqns %%%%%%%%
function output = cAMP_secdegdif(prv)

    k1 = 5;
    k2 = 0.1;
    DcAMP = 50;

    cAMP = prv{1};
    width = size(cAMP,1);
    prv{2} = pdc_calc(prv{2},width);

    cellx = round(prv{2}(:,1));
    celly = round(prv{2}(:,2));
    cell = cAMP*0;

    for i = 1:length(cellx)
            cell(celly(i),cellx(i)) = 1;
    end
    
    try
        output = k1*cell - k2*cAMP + Diffusion(cAMP,DcAMP);
    catch
        [size(cell) size(cAMP)]
    end
end

function output = cellmig(prv)

    f0 = 10;
    a = 1;
    b = 0;
    epsilon = 10^-16;
    sigma = 25;

    cAMP = prv{1};
    width = size(cAMP,1);
    prv{2} = pdc_calc(prv{2},width);

    %%%% gradient of cAMP field %%%%
    cellx = round(prv{2}(:,1));
    celly = round(prv{2}(:,2));
    cx = (circshift(cAMP,[0 -1]) - circshift(cAMP,[0 1]))/2;
    cy = (circshift(cAMP,[-1 0]) - circshift(cAMP,[1 0]))/2;
    % Ncxcy = sqrt(cx.^2 + cy.^2);

    %%%% Lennard-Jones potential %%%%
    celln = length(prv{2}(:,1));
    dx = repmat(prv{2}(:,1),[1 celln]);
    dx = dx' - dx;
    dx(abs(dx) > width/2) = -sign(dx(abs(dx) > width/2)).*(width - abs(dx(abs(dx) > width/2)));
    dy = repmat(prv{2}(:,2),[1 celln]);
    dy = dy' - dy;
    dy(abs(dy) > width/2) = -sign(dy(abs(dy) > width/2)).*(width - abs(dy(abs(dy) > width/2)));
    R = sqrt(dx.^2 + dy.^2);
    LJ = -4*epsilon.*((12*sigma^12)./(R.^13)-(6*sigma^6)./(R.^7));
    LJx = dx./R.*LJ;
    LJy = dy./R.*LJ;
    LJx(1:celln+1:celln^2) = NaN;
    LJy(1:celln+1:celln^2) = NaN;
    LJx = nansum(LJx,2);
    LJy = nansum(LJy,2);
    
    %%%% contact following umeda inouye 2002 %%%%
    velx = repmat(prv{3}(:,1),[1 celln])'./R;
    vely = repmat(prv{3}(:,2),[1 celln])'./R;
    velx = nanmean(velx,2);
    velx(isnan(velx)) = 0;
    vely = nanmean(vely,2);
    vely(isnan(vely)) = 0;
    % vel = sqrt(velx.^2 + vely.^2);

    %%%% contact following original %%%%
    


    %%%% calculation %%%%
    output = prv{2}*0;
    for i = 1:length(cellx)
        motiveFx = cx(celly(i),cellx(i)) + b*velx(i);
        motiveFy = cy(celly(i),cellx(i)) + b*vely(i);
        Norm = sqrt(motiveFx^2 + motiveFy^2);
        output(i,1) = 1/a*(f0*motiveFx/Norm + LJx(i));
        output(i,2) = 1/a*(f0*motiveFy/Norm + LJy(i));
    end
    output(isnan(output)) = 0;
end


%%%%%%%% utility functions %%%%%%%%
function output = Diffusion(input,Dcoeff)

    output = Dcoeff*(circshift(input,[ 1, 0])...
                    +circshift(input,[-1, 0])...
                    +circshift(input,[ 0, 1])...
                    +circshift(input,[ 0,-1])...
                    -4*input);
end

function output = celladd(arg1,varargin)
    varn = length(arg1);
    output = arg1;
    
    for n = 1:length(varargin)
        arg2 = varargin{n};
        for i = 1:varn
            if iscell(arg2)
                output{i} = output{i} + arg2{i};
            else
                output{i} = output{i} + arg2;
            end
        end
    end
end

function output = cellmul(arg1,varargin)
    varn = length(arg1);
    output = arg1;
    
    for n = 1:length(varargin)
        arg2 = varargin{n};
        for i = 1:varn
            if iscell(arg2)
                output{i} = output{i} .* arg2{i};
            else
                output{i} = output{i} * arg2;
            end
        end
    end
end

function output = pdc_calc(arg,width)
    output = arg;
    output(output >= width+0.5) = output(output >= width+0.5) - width;
    output(output < 0.5) = width + output(output < 0.5);
end





