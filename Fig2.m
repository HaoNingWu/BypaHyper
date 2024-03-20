% Bypassing the quadrature exactness assumption of hyperinterpolation on the sphere
% by C. An and H.-N. Wu
% Journal of Complexity (80) 2024, 101789

% MATLAB Codes written by H.-N. Wu in 2022

% Please add the sphere_approx_toolbox_v3.0 onto path

clear
close all

Trial = 10;
ms = [2000:4000:10000,20000,30000,50000:25000:100000,150000:50000:200000];
Ls = [3,6,9,12];

% Validation point set
Xt = get_Xt( );

% Colormap
trad = [46/255, 89/255, 167/255;
    230/255, 0, 18/255;
    250/255 192/255 61/255;
    93/255 163/255 157/255;
    205/255 209/255 113/255;
    240/255 145/255 160/255;
    69/255 70/255 94/255];
figure
rgbplot(trad)
hold on
colormap(trad)
colorbar('Ticks',[])

for func_idx = 2:3
    % function to be approximated
    switch func_idx
        case 1
            func = @(x,y,z) (x+y+z).^2;
            ft = func(Xt(1,:),Xt(2,:),Xt(3,:)); ft = ft';
        case 2
            func = @(x,y,z) (sin(1+abs(x+y+z))).^2+abs(x+y+z);
            ft = func(Xt(1,:),Xt(2,:),Xt(3,:)); ft = ft';
        case 3
            func = @(x,y,z) 0.75*exp(-((9*x-2).^2)/4-((9*y-2).^2)/4-((9*z-2).^2)/4) ...
                +0.75*exp(-((9*x+1).^2)/49-((9*y+1))/10-((9*z+1))/10)...
                +0.5*exp(-((9*x-7).^2)/4-((9*y-3).^2)/4-((9*z-5).^2)/4)...
                -0.2*exp(-((9*x-4).^2)-((9*y-7).^2)-((9*z-5).^2));
            ft = func(Xt(1,:),Xt(2,:),Xt(3,:)); ft = ft';
        case 4
            funtxt_ce = {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'};
            i_f = 3;
            funtxt = funtxt_ce{i_f};
            switch funtxt
                case 'nrWend_k0'
                    rbf_k = 0;
                case 'nrWend_k1'
                    rbf_k = 1;
                case 'nrWend_k2'
                    rbf_k = 2;
                case 'nrWend_k3'
                    rbf_k = 3;
                case 'nrWend_k4'
                    rbf_k = 4;
            end
            switch funtxt
                case {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'}
                    func = @rbf_nr;
            end
            k = rbf_k;
            delta = (3*k+3)*gamma(k+1/2)/(2*gamma(k+1));
            ft = func(Xt',rbf_k);
    end

    for i = 1:max(size(ms))
        m = ms(i);
        for k = 1:Trial
            rvals = 2*rand(m,1)-1;
            elevation = asin(rvals);
            % elevation = -pi/2  + pi*rand(m,1);
            azimuth = 2*pi*rand(m,1);
            radii = ones(m,1);
            [xrand,yrand,zrand] = sph2cart(azimuth,elevation,radii);
            X_k = [xrand,yrand,zrand];

            % function sampling
            switch func_idx
                case {1,2,3}
                    f=func(X_k(:,1),X_k(:,2),X_k(:,3));
                case 4
                    f = func(X_k,rbf_k);
            end

            for ell = 1:max(size(Ls))
                L = Ls(ell);
                fprintf(['(m,n,trial) = (',num2str(m),', ',num2str(L),',',num2str(k),')\n'])

                Y_L = getQ( X_k, L )';
                alpha = 4*pi*Y_L*f/m;

                % approximation polynomial
                L = sqrt(length(alpha))-1;
                Yt_eqp = get_Yt( L, Xt );
                pt_eqphyper = Yt_eqp * alpha;
                errortemp(k,ell) = sqrt(4*pi/50000*sum((ft-pt_eqphyper).^2));
            end
        end

        error(i,:) = sum(errortemp)/Trial;
        errortemp = [];
    end

    % Plotting
    subplot(1,2,func_idx-1)
    semilogy(ms,error(:,1),'marker','square', 'LineWidth',2,...
        'Color',trad(1,:),...
        'MarkerSize',10,...
        'MarkerEdgeColor',trad(1,:),...
        'MarkerFaceColor',trad(1,:))
    hold on
    semilogy(ms,error(:,2),'marker','pentagram', 'LineWidth',2,...
        'Color',trad(2,:)	,...
        'MarkerSize',10,...
        'MarkerEdgeColor',trad(2,:)	,...
        'MarkerFaceColor',trad(2,:))
    semilogy(ms,error(:,3),'marker','o', 'LineWidth',2,...
        'Color',trad(3,:)		,...
        'MarkerSize',10,...
        'MarkerEdgeColor',trad(3,:)	,...
        'MarkerFaceColor',trad(3,:)	)

    semilogy(ms,error(:,4),'marker','d', 'LineWidth',2,...
        'Color',trad(4,:)			,...
        'MarkerSize',10,...
        'MarkerEdgeColor',trad(4,:)		,...
        'MarkerFaceColor',trad(4,:)	)

    if func_idx == 2
        semilogy(ms,4*ms.^(-1/4),'color', trad(5,:), 'LineWidth',2.5)
    else
        semilogy(ms,ms.^(-1/4),'color', trad(5,:), 'LineWidth',2.5)
    end

    set(gca, 'fontsize', 27)
    xlabel('number $m$ of quad. points','interpreter','latex','fontsize',30), ylabel('error','interpreter','latex','fontsize',36)
    box on, grid on,
end

subplot(1,2,1)
title('Convergence of $\mathcal{U}_nf_2$','interpreter','latex','fontsize',36)
set(legend('$n = 3$', '$n=6$', '$n = 9$', '$n=12$','$4m^{-1/4}$'),'interpreter','latex','fontsize',30)

subplot(1,2,2)
title('Convergence of $\mathcal{U}_nf_3$','interpreter','latex','fontsize',36)
set(legend('$n = 3$', '$n=6$', '$n = 9$', '$n=12$','$m^{-1/4}$'),'interpreter','latex','fontsize',30)


