% Bypassing the quadrature exactness assumption of hyperinterpolation on the sphere
% by C. An and H.-N. Wu
% Journal of Complexity (80) 2024, 101789

% MATLAB Codes written by H.-N. Wu in 2022

% Please add the sphere_approx_toolbox_v3.0 onto path

clear
close all

Ls = [6,9,12,15];
ms = [2000:4000:10000,20000,30000,50000:25000:100000,150000:50000:450000];

% Validation point set
Xt = get_Xt( );

% function to be approximated
func = @(x,y,z) (x+y+z).^2;
ft = func(Xt(1,:),Xt(2,:),Xt(3,:)); ft = ft';

for ell = 1:max(size(Ls));
    L = Ls(ell);
    for i = 1:max(size(ms))
        m = ms(i);
        fprintf(['(m,n) = (',num2str(m),', ',num2str(L),')\n'])

        % Equal area points
        X_k = sca_data_cap(2, 1 ,m, 0);
        X_k = X_k';

        % function sampling
        f=func(X_k(:,1),X_k(:,2),X_k(:,3));

        Y_L = getQ( X_k, L )';
        alpha = 4*pi*Y_L*f/m;

        % approximation polynomials
        L = sqrt(length(alpha))-1;
        Yt_eqp = get_Yt( L, Xt );
        pt_eqphyper = Yt_eqp * alpha;
        error(i,ell) = sqrt(4*pi/50000*sum((ft-pt_eqphyper).^2));
    end
end

semilogy(ms,error(:,1),'marker','square', 'LineWidth',2,...
    'Color',[0 0.4470 0.7410],...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0 0.4470 0.7410],...
    'MarkerFaceColor',[0 0.4470 0.7410])
hold on
semilogy(ms,error(:,2),'marker','pentagram', 'LineWidth',2,...
    'Color',[0.8500 0.3250 0.0980]	,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980]	,...
    'MarkerFaceColor',[0.8500 0.3250 0.0980]	)
semilogy(ms,error(:,3),'marker','o', 'LineWidth',2,...
    'Color',[0.9290 0.6940 0.1250]		,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.9290 0.6940 0.1250]	,...
    'MarkerFaceColor',[0.9290 0.6940 0.1250]	)
semilogy(ms,error(:,4),'marker','d', 'LineWidth',2,...
    'Color',[0.4660 0.6740 0.1880]		,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.4660 0.6740 0.1880]			,...
    'MarkerFaceColor',[0.4660 0.6740 0.1880]		)

semilogy(ms,ms.^(-1/2),'-k', 'LineWidth',2.5)
semilogy(ms,ms.^(-1),'--k', 'LineWidth',2.5)

set(gca, 'fontsize', 27)
xlabel('number $m$','interpreter','latex','fontsize',33), ylabel('error','interpreter','latex','fontsize',42)
set(legend( '$n=6$','$n=9$', '$n = 12$', '$n=15$','$m^{-1/2}$','$m^{-1}$'),'interpreter','latex','fontsize',30)
box on, grid on