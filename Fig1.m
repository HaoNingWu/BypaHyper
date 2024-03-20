% Bypassing the quadrature exactness assumption of hyperinterpolation on the sphere
% by C. An and H.-N. Wu
% Journal of Complexity (80) 2024, 101789

% MATLAB Codes written by H.-N. Wu in 2022

% Please add the sphere_approx_toolbox_v3.0 onto path

clear
close all

Trial = 10;
ms = [500,1000,2000:4000:10000,20000,30000,50000:25000:100000,150000:50000:250000];
Ls = [6,9,12,15];

% Validation point set
Xt = get_Xt( );

func_idx = 1;
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

% Testing started
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
    errorub(i,:) = max(errortemp);
    errorlb(i,:) = min(errortemp);
    error(i,:) = sum(errortemp)/Trial;
    errortemp = []; etalbtemp = [];
end

% Plotting error curves
semilogy(ms,error(:,1),'marker','square', 'LineWidth',2,...
    'Color','b',...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b')
hold on
semilogy(ms,error(:,2),'marker','pentagram', 'LineWidth',2,...
    'Color','r',...
    'MarkerSize',10,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
semilogy(ms,error(:,3),'marker','o', 'LineWidth',2,...
    'Color','y',...
    'MarkerSize',10,...
    'MarkerEdgeColor','y',...
    'MarkerFaceColor','y')
semilogy(ms,error(:,4),'marker','d', 'LineWidth',2,...
    'Color','g',...
    'MarkerSize',10,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','g')
semilogy(ms,ms.^(-1/4),'color', 'k', 'LineWidth',2.5)

% Plotting error bands
Gd = [ms flip(ms)];
Bd1 = [errorub(:,1)' flip(errorlb(:,1)')];
h = fill(Gd,Bd1,'b','facealpha',0.1);
Bd2 = [errorub(:,2)' flip(errorlb(:,2)')];
h = fill(Gd,Bd2,'r','facealpha',0.1);
Bd3 = [errorub(:,3)' flip(errorlb(:,3)')];
h = fill(Gd,Bd3,'y','facealpha',0.1);
Bd4 = [errorub(:,4)' flip(errorlb(:,4)')];
h = fill(Gd,Bd4,'g','facealpha',0.1);

set(gca, 'fontsize', 27)
set(legend('$n = 6$', '$n=9$', '$n = 12$', '$n=15$','$m^{-1/4}$'),'interpreter','latex','fontsize',30)
title('Convergence of $\mathcal{U}_nf_1$','interpreter','latex','fontsize',32)
xlabel('number $m$ of quad. points','interpreter','latex','fontsize',30), ylabel('error','interpreter','latex','fontsize',36)
box on, grid on,
