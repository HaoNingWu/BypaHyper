% Bypassing the quadrature exactness assumption of hyperinterpolation on the sphere
% by C. An and H.-N. Wu
% Journal of Complexity (80) 2024, 101789

% MATLAB Codes written by H.-N. Wu in 2022

% Please add the sphere_approx_toolbox_v3.0 onto path

clear
close all

L = 5;
point_idx = 1;
ms = [4           16          36          64         100         121         441         961        1681        2601        3721        5041        6561  8281       10201];
Ks = [0,1,2,3,4];

% Validation point set
Xt = get_Xt( );

for ell = 1:max(size(Ks));
    i_f = Ks(ell);

    for i = 1:max(size(ms))
        m = ms(i);
        fprintf(['(pointset_idx,m,sigma) = (',num2str(point_idx),', ',num2str(m),', ',num2str(i_f),')\n'])

        % function to be approximated
        funtxt_ce = {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'};
        funtxt = funtxt_ce{i_f+1};
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

        switch point_idx
            case 1 % Equal area points
                X_k = sca_data_cap(2, 1 ,m, 0);
                X_k = X_k';
            case 2 % Minimal energy points
                X_k = loadME( t, (t+1)^2 );
            case 3
                X_k = loadMD( t, (t+1)^2 );
            case 4
                X_k = loadStd(t,(t+1)^2);
        end

        % function sampling
        f = func(X_k,rbf_k);
        Y_L = getQ( X_k, L )';
        alpha = 4*pi*Y_L*f/m;

        % approximation polynomials
        L = sqrt(length(alpha))-1;
        Yt_eqp = get_Yt( L, Xt );
        pt_eqphyper = Yt_eqp * alpha;
        error(i,ell) = sqrt(4*pi/50000*sum((ft-pt_eqphyper).^2));
    end
end


axes('position',[0.07 0.12 0.18 0.8]),
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
semilogy(ms,error(:,5),'marker','^', 'LineWidth',2,...
    'Color',[0.3010 0.7450 0.9330]			,...
    'MarkerSize',10,...
    'MarkerEdgeColor',[0.3010 0.7450 0.9330]				,...
    'MarkerFaceColor',[0.3010 0.7450 0.9330]			)

set(gca, 'fontsize', 27)
xlabel('number $m$','interpreter','latex','fontsize',33), ylabel('error','interpreter','latex','fontsize',42)
set(legend( '$\sigma=0$','$\sigma=1$', '$\sigma=2$', '$\sigma=3$', '$\sigma=4$'),'interpreter','latex','fontsize',30)
title('Equal area points','interpreter','latex','fontsize',36)
box on, grid on,


ts = [1:2:9,10:10:100];
Ks = [0,1,2,3,4];
for point_idx = 2:4
    for ell = 1:max(size(Ks));
        i_f = Ks(ell);
        for i = 1:max(size(ts))
            t = ts(i);
            fprintf(['(pointset_idx,t,sigma) = (',num2str(point_idx),', ',num2str(t),', ',num2str(i_f),')\n'])

            % function to be approximated
            funtxt_ce = {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'};
            funtxt = funtxt_ce{i_f+1};
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

            switch point_idx
                case 1 % Equal area points
                    X_k = sca_data_cap(2, 1 ,m, 0);
                    X_k = X_k';
                case 2 % Minimal energy points
                    X_k = loadME( t, (t+1)^2 );
                case 3
                    X_k = loadMD( t, (t+1)^2 );
                case 4
                    X_k = loadStd(t,(t+1)^2);
            end

            [m,n] = size(X_k);
            ms(i) = m;

            %function sampling
            f = func(X_k,rbf_k);
            Y_L = getQ( X_k, L )';
            alpha = 4*pi*Y_L*f/m;

            % approximation polynomials
            L = sqrt(length(alpha))-1;
            Yt_eqp = get_Yt( L, Xt );
            pt_eqphyper = Yt_eqp * alpha;
            error(i,ell) = sqrt(4*pi/50000*sum((ft-pt_eqphyper).^2));
        end
    end

    axes('position',[0.07+0.23*(point_idx-1) 0.12 0.18 0.8]),
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
    semilogy(ms,error(:,5),'marker','^', 'LineWidth',2,...
        'Color',[0.3010 0.7450 0.9330]		,...
        'MarkerSize',10,...
        'MarkerEdgeColor',[0.3010 0.7450 0.9330]			,...
        'MarkerFaceColor',[0.3010 0.7450 0.9330]	)
    set(gca, 'fontsize', 27)
    xlabel('number $m$','interpreter','latex','fontsize',33),

    set(legend( '$\sigma=0$','$\sigma=1$', '$\sigma=2$', '$\sigma=3$', '$\sigma=4$'),'interpreter','latex','fontsize',30)
    box on, grid on,

    if point_idx == 2
        title('Coulomb energy points','interpreter','latex','fontsize',36)
    elseif point_idx == 3
        title('Fekete points','interpreter','latex','fontsize',36)
    else
        title('Spherical $t$-designs','interpreter','latex','fontsize',36)
    end

end
