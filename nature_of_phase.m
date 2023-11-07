figure;
axis off;
text(0, 1,'The Nature of Phase','Units','normalized', 'fontsize', 20, 'Interpreter','latex')
text(0, 0.9,'$X(n,k)$: $n$ frame index, $k$ bin index with $0$ for DC bin; $\angle$: takes phase','Units','normalized', 'fontsize', 15, 'Interpreter','latex')
text(0, 0.8,'(1): $\angle  E\left[ { X(n+1, k) }{X^*(n, k) } \right] = 2\pi k\frac{\rm block\_size}{{\rm fft\_size}}$','Units','normalized', 'fontsize', 15, 'Interpreter','latex')
text(0, 0.7,'(2): $\angle E\left[ X(n,k+1)X^*(n,k) \right] = 2\pi \frac{{ \rm group\_delay\_of\_}k{\rm th\_analysis\_filter}}{\rm fft\_size}  $','Units','normalized', 'fontsize', 15, 'Interpreter','latex')
text(0, 0.6,'(3): Nonexistence of an analytic structural phase bias function $\phi(n,k)$: ','Units','normalized', 'fontsize', 15, 'Interpreter','latex')
text(0, 0.5,'$\qquad$ Let $\frac{\partial \phi(n,k)}{\partial } $ takes finite difference wrt $n$ or $k$','Units','normalized', 'fontsize', 15, 'Interpreter','latex')
text(0, 0.4,'$\qquad$ Eq. (1) sugguests $\frac{\partial{ \phi_1(n,k)}}{\partial n}=2\pi k\frac{\rm blockq\_size}{{\rm fft\_size}}$','Units','normalized', 'fontsize', 15, 'Interpreter','latex')
text(0, 0.3,'$\qquad$ Eq. (2) sugguests $\frac{\partial{ \phi_2(n,k)}}{\partial k}=2\pi \frac{ \rm group\_delay }{\rm fft\_size}$','Units','normalized', 'fontsize', 15, 'Interpreter','latex')
text(0, 0.2,'$\qquad$ Generally, $\frac{\partial}{\partial k} \frac{\partial{ \phi_1(n,k)}}{\partial n} \ne \frac{\partial}{\partial n} \frac{\partial{ \phi_2(n,k)}}{\partial k} $','Units','normalized', 'fontsize', 15, 'Interpreter','latex')
text(0, 0.1,'$\qquad$ $\Rightarrow$ inconsistent PDE!','Units','normalized', 'fontsize', 15, 'Interpreter','latex')