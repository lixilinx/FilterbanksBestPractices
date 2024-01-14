figure;
axis off;
text(0, 1,'Designs with different smallest momentums','Units','normalized', 'fontsize', 15, 'Interpreter','latex')

text(0, 11/12,'min 0th order momentum (default)','Units','normalized', 'fontsize', 12, 'Interpreter','latex')
text(0, 10/12,'$\int _{\omega _c}^{\pi } \left|\sum _{t=0}^{L - 1} h_p e^{-\jmath \omega p} \right|^2 d\omega = \bf {h}^{T} \,\bf {\Pi }_0 \,\bf h $ with $({\bf \Pi }_0)_{p,\, q}$ be','Units','normalized', 'fontsize', 12, 'Interpreter','latex')
text(0, 9/12,'$\pi-\omega_c$ if $p=q$ else $-\frac{\sin[(p-q)\omega_c]}{p-q}$','Units','normalized', 'fontsize', 12, 'Interpreter','latex')

text(0, 7/12,'min 1st order momentum','Units','normalized', 'fontsize', 12, 'Interpreter','latex')
text(0, 6/12,'$\int _{\omega _c}^{\pi } \omega \left|\sum _{t=0}^{L - 1} h_p e^{-\jmath \omega p} \right|^2 d\omega = \bf {h}^{T} \,\bf {\Pi }_1 \,\bf h $ with $({\bf \Pi }_1)_{p,\, q}$ be','Units','normalized', 'fontsize', 12, 'Interpreter','latex')
text(0, 5/12,'$\frac{\pi^2-\omega_c^2}{2}$ if $p=q$ else $-\frac{\cos[(p-q)\omega_c] + (p-q)\omega_c \sin[(p-q)\omega_c] - (-1)^{p-q}}{(p-q)^2}$','Units','normalized', 'fontsize', 12, 'Interpreter','latex')

text(0, 3/12,'min 2nd order momentum','Units','normalized', 'fontsize', 12, 'Interpreter','latex')
text(0, 2/12,'$\int _{\omega _c}^{\pi } \omega^2 \left|\sum _{t=0}^{L - 1} h_p e^{-\jmath \omega p} \right|^2 d\omega = \bf {h}^{T} \,\bf {\Pi }_2 \,\bf h $ with $({\bf \Pi }_2)_{p,\, q}$ be','Units','normalized', 'fontsize', 12, 'Interpreter','latex')
text(0, 1/12,'$\frac{\pi^3-\omega_c^3}{3}$ if $p=q$ else $-\frac{[(p-q)^2\omega_c^2 - 2]\sin[(p-q)\omega_c] + 2(p-q)\omega_c\cos[(p-q)\omega_c] - 2\pi(p-q) (-1)^{p-q}}{(p-q)^3}$','Units','normalized', 'fontsize', 12, 'Interpreter','latex')