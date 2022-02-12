"""
Convert the elastic thickness into
a heat flow given the input rheology.
Convert a temperature profile into
a yield strength envelope and elastic
thickness.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def Curv_Moment(
    d_sig_tmp0,
    d_sig_tmp1,
    K_curv,
    young,
    poisson,
    depth,
    step_depth,
    sig_y,
    HF=0,
    neutralfib=None,
    plot_YSE=False,
    decoupling=False,
    show=True,
    ax=None,
    figsize=mpl.rcParams["figure.figsize"],
):
    """
    Determine the bending moment given the input yield
    strength envelope and curvature.

    Returns
    -------
    Mx : float
       The bending moment of the plate
    d_sig_tab1 : array size(depth)
       Integrated part of the yield strength envelope in compression (Pa)
    d_sig_tab2 : array size(depth)
       Integrated part of the yield strength envelope in tension (Pa)
    sig_ela : array size(depth)
       Integrated elastic part of the yield strength envelope (Pa)
    neutralfib : float
        Neutral fiber index to be used for the next iteration

    Parameters
    ----------
    d_sig_tmp0 : float
       Stress curve, tension (Pa).
    d_sig_tmp1 : float
       Stress curve, compression (Pa).
    K_curv : float
       Plate curvature (m-1).
    young : float
       Young's modulus (Pa).
    poisson : float
       Poisson's ratio.
    depth : array, float
       Input depth for the integration (m).
    step_depth : float
       Iteration depth step for the integration.
    sig_y : float
       Bounding stress (Pa).
    HF : float, optional, default = 0
       The input mantle heat flow (W m-2).
    neutralfib : int, optional, default = None
        Neutral fiber index of the previous interation.
    plot_YSE : boolean, optional, default = False
       If true, plot the yield strength envelope
    decoupling : boolean, optional, default = False
       If true, the rheology is decoupled and we plot 2 sets
       of curvatures on ax if plot_YSE.
    show : boolean, optional, default = True
       If true, show the plot once the first decoupled part is plotted on
       ax if plot_YSE.
    ax : object, optional, default = None
       Axis for plotting.
    figsize : tuple, optional, default = mpl.rcParams['figure.figsize']
       Size of the output figure.
    """

    first = False
    if neutralfib is None:
        for i in depth:
            # Find neutral fiber
            sig_ela = (-young * K_curv * (depth - i)) / (1.0 - poisson**2)
            d_sig_tab1 = np.minimum(sig_ela[:], d_sig_tmp0[:])
            d_sig_tab2 = np.maximum(sig_ela[:], d_sig_tmp1[:])
            if K_curv >= 0:
                # elastic core limit & neutral fiber
                bound2 = np.min(np.where(sig_ela == 0))
                d_sig_tab1[bound2:] = 0.0
                d_sig_tab2[:bound2] = 0.0

                NetAxial = np.trapz(d_sig_tab1[:bound2], dx=step_depth) + np.trapz(
                    d_sig_tab2[bound2:], dx=step_depth
                )
            else:
                # elastic core limit & neutral fiber
                bound2 = np.min(np.where(sig_ela == 0))
                d_sig_tab1[:bound2] = 0.0
                d_sig_tab2[bound2:] = 0.0

                NetAxial = np.trapz(d_sig_tab2[:bound2], dx=step_depth) + np.trapz(
                    d_sig_tab1[bound2:], dx=step_depth
                )

            if not first and NetAxial != 0:
                tempAx = NetAxial
                first = True
            if first and tempAx < 0 and NetAxial > 0:
                neutralfib = i
                break
            elif first and tempAx > 0 and NetAxial < 0:
                neutralfib = i
                break
    else:
        sig_ela = (-young * K_curv * (depth - neutralfib)) / (1.0 - poisson**2)
        d_sig_tab1 = np.minimum(sig_ela[:], d_sig_tmp0[:])
        d_sig_tab2 = np.maximum(sig_ela[:], d_sig_tmp1[:])
        if K_curv >= 0:
            # elastic core limit & neutral fiber
            bound2 = np.min(np.where(sig_ela == 0))
            d_sig_tab1[bound2:] = 0.0
            d_sig_tab2[:bound2] = 0.0
        else:
            # elastic core limit & neutral fiber
            bound2 = np.min(np.where(sig_ela == 0))
            d_sig_tab1[:bound2] = 0.0
            d_sig_tab2[bound2:] = 0.0

    d_sig_tab3 = d_sig_tab1 * (depth - neutralfib / 2.0)
    d_sig_tab4 = d_sig_tab2 * (depth - neutralfib / 2.0)

    Mx = np.trapz(d_sig_tab3, dx=step_depth) + np.trapz(d_sig_tab4, dx=step_depth)

    if plot_YSE:  # Show YSE
        depth_km = depth / 1e3
        if not decoupling:
            f, ax = plt.subplots(1, 1, figsize=figsize)
        ax.plot(d_sig_tmp0, depth_km, "purple")
        ax.plot(d_sig_tmp1, depth_km, "k")
        ax.plot(d_sig_tab1, depth_km, color="orange", label="Integrated tension")
        ax.plot(d_sig_tab2, depth_km, "b", label="Integrated compression")
        ax.plot(sig_ela, depth_km, "--r", label="Curvature")
        if show:
            ax.set_xlim(
                1.5 * np.min([d_sig_tmp1.min(), d_sig_tab2.min()]),
                1.5 * np.max([d_sig_tmp0.max(), d_sig_tmp1.max()]),
            )
            ax.axvline(sig_y, label="Bounding stress")
            ax.axvline(-sig_y)
            ax.set_ylabel("Depth (km)")
            ax.set_xlabel("Sress (Pa)")
            ax.legend(loc="lower left")
            ax.invert_yaxis()
            if HF != 0:
                ax.set_title("Mantle heat flow %s (mW m-2)" % (HF * 1e3))
            plt.show()

    return Mx, d_sig_tab1, d_sig_tab2, sig_ela, neutralfib


def Conversion_Te_HF(
    Q_crust,
    A_crust,
    n_crust,
    Q_mantle,
    A_mantle,
    n_mantle,
    g_surf,
    rhobar,
    R_mpr,
    Te,
    Tc,
    k_crust,
    k_mantle,
    K_curv,
    H_c,
    rhoc,
    rhom,
    eps,
    young,
    poisson,
    sig_y,
    Ts,
    R_gas,
    HF_min=None,
    HF_max=None,
    HF_step=None,
    step_depth=None,
    max_depth=None,
    quiet=True,
    plot=False,
    plot_YSE=False,
    figsize=mpl.rcParams["figure.figsize"],
):
    """
    Determine the surface, crustal, and mantle heat flows, mechanical thickness,
    and thermal gradients from input rheology and elastic parameters.

    Returns
    -------
    F_c : float
       The crustal heat flow (mW m-2).
    F_m : float
       The mantle heat flow (mW m-2).
    F_s : float
       The surface heat flow (mW m-2).
    dT_c : float
       The crust temperature gradient (K km-1).
    dT_m : float
       The mantle temperature gradient (K km-1).
    dT_s : float
       The surface temperature gradient (K km-1).
    Tm : Float
       Mechanical thickness of the lithosphere (m).
    d_sig_tmp1_best : array, size(max_depth / step_depth)
       The best-fitting YSE compression (Pa).
    d_sig_tmp0_best : array, size(max_depth / step_depth)
       The best-fitting YSE tension (Pa).
    ax1 : matplotlib axes object
        If plot is True – Axis where the YSE is plotted.
    ax2 : matplotlib axes object
        If plot is True – Axis where the temperature profile is plotted.

    Parameters
    ----------
    Q_crust : float
       Activation energy for the crust (J mol-1).
    A_crust : float
       Empirical constant for the crust (Pa-n s-1).
    n_crust : float
       Rheological constant for the crust.
    Q_mantle : float
       Activation energy for the mantle (J mol-1).
    A_mantle : float
       Empirical constant for the mantle (Pa-n s-1).
    n_mantle : float
       Rheological constant for the mantle.
    g_surf : float
       Gravitational attraction of the surface (m s-2).
    rhobar : float
       Mean density of the body (kg m-3).
    R_mpr : float
       Mean radius of the body (m).
    Te : float
       Elastic thickness (m).
    Tc : float
       Crustal thickness (m).
    k_crust : float
       Thermal conductivity in the crust (W m−1 K−1).
    k_mantle : float
       Thermal conductivity in the mantle (W m−1 K−1).
    K_curv : float
       Plate curvature (m-1).
    H_c : float
       Average volumetric crustal heat production (W kg-1).
    rhoc : float
       Crustal density (kg m-3).
    rhom : float
       Mantle density (kg m-3).
    eps : float
       Strain rate (s-1).
    young : float
       Young's modulus (Pa).
    poisson : float
       Poisson's ratio.
    sig_y : float
       Bounding stress (Pa).
    Ts : float
       Surface temperature (K).
    R_gas : float
       Gas constant (J mol-1 K-1).
    HF_min : float, optional, default = None
       Minimum heat flow to test (mW m-2), if None
       we guess it assuming mechanical thickness
       equals elastic thickness (zero curvature).
    HF_max : float, optional, default = None
       Maximum heat flow to test (mW m-2), if None
       we guess it assuming mechanical thickness
       equals elastic thickness (zero curvature).
    HF_step : float, optional, default = 1
       Heat flow step to test (mW m-2)
    step_depth : integer, optional, default = Te / 1000.
       Iteration depth step for the integration.
    max_depth : integer, optional, default = Te * 3
       Maximum depth for the integration.
    quiet : boolean, optional, default = True
       if True, print various outputs.
    plot : boolean, optional, default = False
       if True, plot the best-fit yield strength envelope
       and temperature profile.
    plot_YSE : boolean, optional, default = False
       if True, plot the yield strength envelope for each
       tested heat flows.
    figsize : tuple, optional, default = mpl.rcParams['figure.figsize']
       Size of the output figure.
    """

    if not quiet:
        print("Elastic thickness (km): %i" % (Te / 1e3))
        print("Crustal thickness (km): %i" % (Tc / 1e3))
        print("mantle and crustal densities (kg m-3): %i, %i" % (rhom, rhoc))
        print(
            "mantle and crustal thermal conductivities (W m−1 K−1): %s, %s"
            % (k_mantle, k_crust)
        )
        print("Bounding stress (MPa): %i" % (sig_y / 1e6))
        print("Crustal heat production (W kg-1): %s" % (H_c))
        print("Curvature (m-1): %s" % (K_curv))
        print("Strain rate (s-1): %s" % (eps))

    # First order heat flow calculation
    # Valid when the curvature is low, about 1e-9 (see McNutt 1984)
    # Doesn't account for a potential decoupling between the crust and mantle,
    # which can vary the elastic thickness by about 30%, hence the derived heat flow
    k_avg = (k_mantle * (Te - Tc) + k_crust * Tc) / Te  # Linear avg
    if Te > Tc:  # Mantle rheology
        Tmax = (Q_mantle / R_gas) / (np.log(A_mantle * (sig_y**n_mantle) / eps))
        F = 1e3 * k_avg * (Tmax - Ts) / Te
    else:  # Crustal rheology
        Tmax = (Q_crust / R_gas) / (np.log(A_crust * (sig_y**n_crust) / eps))
        F = 1e3 * k_crust * (Tmax - Ts) / Te

    if not quiet:
        print("Threshold T° at ductile failure ", Tmax)
        print(
            "First order, no crustal heat, surface (= lithospheric) heat flow %.2f mW m-2"
            % (F)
        )
        print(
            "Valid when curvature is low (< 1e-9), but doesn't account for crust/mantle decoupling"
        )

    # Initial guesses on the heat flows to be explored
    if HF_min is None:
        HF_min = int(0.3 * F)
    if HF_max is None:
        HF_max = int(1.5 * F)
    if HF_step is None:
        HF_step = 1
    HF_arr = range(HF_min, HF_max, HF_step)

    if max_depth is None:
        max_depth = np.max(
            [int(Te * 3), int(Tc + 2 * Te)]
        )  # max depth for the integration (m)
    if step_depth is None:
        step_depth = int(
            Te / 1000
        )  # step of the integration, determines the accuracy of the numerical integration
    # but slows down the computation.
    depths_arr = np.array(range(1, int(max_depth), step_depth))

    d_sig_tab = np.zeros((3, np.size(depths_arr)))
    T_z_tab = np.zeros_like(depths_arr)

    # Integration of the elastic strength of the lithosphere
    d_sig_tab[2, :] = -(K_curv * young) / float(1.0 - poisson**2)  # Pa
    a = -1
    for j in depths_arr:
        a = a + 1
        d_sig_tab[2, a] = (
            -d_sig_tab[2, a] * (j - float(Te) / 2.0) ** 2
        )  # To integrate the elastic strength
        if j > Te:
            bound4 = a
            break
    M_el = np.trapz(
        d_sig_tab[2, :bound4], dx=step_depth
    )  # The analytic formula for M_el is not used here.
    # as we use trapz for M_real, we use it to compute M_el
    # but in practice M_el integrated and the analytic formula are close
    M_el_analyt = K_curv * young * Te**3 / (12.0 * float(1.0 - poisson**2))
    if not quiet:
        print(
            "Numerical integration of M_el vs exact integration M_el",
            "%.4e" % (M_el),
            "%.4e" % (M_el_analyt),
        )

    if np.abs(M_el - M_el_analyt) / M_el_analyt > 5e2:
        raise ValueError(
            "The difference between the integrated elastic moment"
            + " and analytical is too large (%s)"
            % (100 * np.abs(M_el - M_el_analyt) / M_el_analyt)
            + "Increase either step_depth or max_depth, which are here set to %i and %i km"
            % (step_depth / 1e3, max_depth / 1e3)
        )

    misfit_temp = 1e50
    decoupling_potential = False
    decoupling = False
    prints_weak = False
    prints_decoup = False
    prints_decoup_m = False

    d_sig_tab0B, d_sig_tab1B = Brittle_Strength(
        Tc, R_mpr, g_surf, rhom, rhoc, rhobar, depths_arr
    )

    for HF in HF_arr:  # HF : mw/m2
        d_sig_tab[0, :] = d_sig_tab0B
        d_sig_tab[1, :] = d_sig_tab1B
        duct_fail_c = False
        duct_fail_m = False
        F = float(HF / 1e3)
        for iter, z in enumerate(depths_arr):
            # Ductile
            # Temperature
            if z <= Tc:
                Fs = F + rhoc * H_c * z
                T_z = Ts + Fs * z / k_crust - rhoc * H_c * z**2 / (2 * k_crust)

                # Onset of ductile deformation
                ductile_sig = (eps / A_crust) ** (1.0 / float(n_crust)) * np.exp(
                    Q_crust / (n_crust * R_gas * T_z)
                )
                d_sig_tab[0, iter] = min(ductile_sig, d_sig_tab[0, iter])
                d_sig_tab[1, iter] = -min(ductile_sig, abs(d_sig_tab[1, iter]))

                # Crustal ductile failure
                if (
                    d_sig_tab[0, iter] == ductile_sig
                    or d_sig_tab[1, iter] == -ductile_sig
                ) and not duct_fail_c:
                    duct_fail_c = True
                    i_duct_c = iter
                    i_crust_thick = np.argmin((depths_arr - Tc) ** 2)
            else:
                T_z_c = Ts + Fs * Tc / k_crust - rhoc * H_c * Tc**2 / (2 * k_crust)
                T_z = T_z_c + F / k_mantle * (z - Tc)

                # Onset of ductile deformation
                ductile_sig = (eps / A_mantle) ** (1.0 / float(n_mantle)) * np.exp(
                    Q_mantle / (n_mantle * R_gas * T_z)
                )
                d_sig_tab[0, iter] = min(ductile_sig, d_sig_tab[0, iter])
                d_sig_tab[1, iter] = -min(ductile_sig, abs(d_sig_tab[1, iter]))

                # Mantle ductile failure
                if (
                    d_sig_tab[0, iter] == ductile_sig
                    or d_sig_tab[1, iter] == -ductile_sig
                ) and not duct_fail_m:
                    duct_fail_m = True
                    i_duct_m = iter

            T_z_tab[iter] = T_z

        # Bounding stress limit
        if duct_fail_m:
            depth_m0 = np.where(abs(d_sig_tab[0, i_duct_m:]) <= sig_y)
            depth_m1 = np.where(abs(d_sig_tab[1, i_duct_m:]) <= sig_y)
            d_sig_tab[0, i_duct_m:][depth_m0] = 0.0
            d_sig_tab[1, i_duct_m:][depth_m1] = 0.0
        if duct_fail_c:
            depth_c0 = np.where(abs(d_sig_tab[0, i_duct_c:i_crust_thick]) <= sig_y)
            depth_c1 = np.where(abs(d_sig_tab[1, i_duct_c:i_crust_thick]) <= sig_y)
            d_sig_tab[0, i_duct_c:i_crust_thick][depth_c0] = 0.0
            d_sig_tab[1, i_duct_c:i_crust_thick][depth_c1] = 0.0

        if not decoupling_potential:
            M_real, d_sig_tab1, d_sig_tab2, sig_ela, neutralfib_nd = Curv_Moment(
                d_sig_tab[0, :],
                d_sig_tab[1, :],
                K_curv,
                young,
                poisson,
                depths_arr,
                step_depth,
                sig_y,
                show=False,
                decoupling=False,
            )
            ### Determine whether the rheology is decoupled.
            if duct_fail_c:  # Ductile failure in the crust required for decoupling
                if (
                    d_sig_tab[0, i_crust_thick:] > 0
                ).any() and (  # Do we have strength below the crust
                    (
                        sig_ela[i_duct_c:][depth_m0] > d_sig_tab[0, i_duct_c:][depth_m0]
                    ).any()  # Is the compressional elastic strength larger than the ductile failure in the crust
                    or (
                        sig_ela[i_duct_c:][depth_m1] < d_sig_tab[1, i_duct_c:][depth_m1]
                    ).any()  # Is the extensional elastic strength larger than the ductile failure in the crust
                ):
                    decoupling_potential = True

            if duct_fail_c and not quiet and not prints_weak:
                print("## Apparition of a weak ductile upper layer ##")
                prints_weak = True

        if decoupling_potential:
            d_sig_tab_c0 = d_sig_tab[0, :].copy()
            d_sig_tab_c1 = d_sig_tab[1, :].copy()
            d_sig_tab_m0 = d_sig_tab[0, :].copy()
            d_sig_tab_m1 = d_sig_tab[1, :].copy()

            i_duct_base_c0 = np.argmin(abs(d_sig_tab_c0[i_duct_c:]))
            i_duct_base_c1 = np.argmin(abs(d_sig_tab_c1[i_duct_c:]))
            d_sig_tab_c0[i_duct_c:][i_duct_base_c0:] = 0.0
            d_sig_tab_c1[i_duct_c:][i_duct_base_c1:] = 0.0

            d_sig_tab_m0[: i_duct_base_c0 + i_duct_c] = 0.0
            d_sig_tab_m1[: i_duct_base_c1 + i_duct_c] = 0.0
            sum_d0_m = np.sum(d_sig_tab_m0) != 0
            sum_d1_m = np.sum(d_sig_tab_m1) != 0

            if plot_YSE:
                f, ax = plt.subplots(1, 1, figsize=figsize)
            else:
                ax = None

            args_M_real = (K_curv, young, poisson, depths_arr, step_depth, sig_y)
            args_M_real_d = dict(
                plot_YSE=plot_YSE, HF=F, decoupling=decoupling_potential, ax=ax
            )

            if sum_d0_m and sum_d1_m:
                M_real_m, d_sig_tab1_m, d_sig_tab2_m, sig_ela_m, _ = Curv_Moment(
                    d_sig_tab_m0,
                    d_sig_tab_m1,
                    *args_M_real,
                    **args_M_real_d,
                    show=False,
                )
                decoupling = True
                if not prints_decoup:
                    prints_decoup = True
                    if not quiet:
                        print(
                            "## Beginning of decoupling between the crust and mantle ##"
                        )
            else:
                M_real_m = 0
                decoupling = False
                decoupling_potential = False
                if not quiet and not prints_decoup_m and prints_decoup:
                    print("## Mantle layer has no more strength ##")
                    prints_decoup_m = True

            M_real_c, d_sig_tab1_c, d_sig_tab2_c, sig_ela_c, neutralfib = Curv_Moment(
                d_sig_tab_c0,
                d_sig_tab_c1,
                *args_M_real,
                **args_M_real_d,
                neutralfib=neutralfib_nd if not decoupling else None,
                show=plot_YSE,
            )
            plotted = True
            M_real = M_real_c + M_real_m
        else:
            plotted = False

        if not decoupling:
            M_real, d_sig_tab1, d_sig_tab2, sig_ela, _ = Curv_Moment(
                d_sig_tab[0, :],
                d_sig_tab[1, :],
                K_curv,
                young,
                poisson,
                depths_arr,
                step_depth,
                sig_y,
                neutralfib=neutralfib_nd,
                plot_YSE=not plotted and plot_YSE,
                HF=F,
            )

        misfit = abs(1.0 - abs(M_real / M_el))
        if misfit <= misfit_temp:
            if not quiet:
                print(
                    " YSE misfit %.2f, M_real/M_el %.2f, crust, mantle, and surface heat flows (mW m-2) %.2f, %.2f, %.2f "
                    % (misfit, abs(M_real / M_el), (Fs - F) * 1e3, F * 1e3, Fs * 1e3)
                )
            misfit_temp = misfit
            F_m_best = F * 1e3
            F_s_best = Fs * 1e3
            d_sig_tmp0_best = d_sig_tab[0, :].copy()
            d_sig_tmp1_best = d_sig_tab[1, :].copy()
            if decoupling:
                sig_ela_best_c = sig_ela_c
                sig_ela_best_m = sig_ela_m
                decoupling_better = True
                i_duct_m_best = i_duct_m
                i_duct_c_best = i_duct_c
                depth_m0_best = depth_m0
                d_sig_tab1_best = d_sig_tab1_c + d_sig_tab1_m
                d_sig_tab2_best = d_sig_tab2_c + d_sig_tab2_m
            else:
                sig_ela_best = sig_ela
                decoupling_better = False
                d_sig_tab1_best = d_sig_tab1
                d_sig_tab2_best = d_sig_tab2
            T_z_best = T_z_tab
        else:
            if not quiet:
                print(
                    " YSE misfit is worse %.2f, M_real/M_el %.2f, crust, mantle, and surface heat flows (mW m-2) %.2f, %.2f, %.2f "
                    % (misfit, abs(M_real / M_el), (Fs - F) * 1e3, F * 1e3, Fs * 1e3)
                )

    # Mechanical thickness of the lithosphere
    find_Tm = np.where(np.logical_and(abs(d_sig_tmp0_best) <= sig_y, depths_arr > Te))
    if np.size(find_Tm) == 0:
        raise ValueError(
            " Mechanical thickness larger than maximum input depth\n"
            + "User should increase the max_depth parameter, value is: %i km"
            % (max_depth / 1e3)
        )
    else:
        if decoupling_better:
            Tm1 = depths_arr[find_Tm][0]
            Tm2 = (
                depths_arr[i_duct_m_best:][depth_m0_best][0] - depths_arr[i_duct_m_best]
            )
            Tm = (Tm1**3 + Tm2**3) ** (1 / 3)  # Burov & Diament 1995
        else:
            Tm = depths_arr[find_Tm][0]
    F_c_best = F_s_best - F_m_best

    if not quiet:
        print(
            "Best-fit crust, mantle, and surface heat flows (mW m-2) %.2f, %.2f, %.2f"
            % (F_c_best, F_m_best, F_s_best)
        )
        print("Mechanical thickness of the lithosphere (km) %i" % (Tm / 1e3))

    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        depth_km = depths_arr / 1e3
        ax1.plot(d_sig_tmp0_best, depth_km, "purple")
        ax1.plot(d_sig_tmp1_best, depth_km, "k")
        ax1.plot(d_sig_tab1_best, depth_km, color="orange", label="Integrated tension")
        ax1.plot(d_sig_tab2_best, depth_km, "b", label="Integrated compression")
        if decoupling_better:
            ax1.plot(
                sig_ela_best_c[:i_duct_m_best],
                depth_km[:i_duct_m_best],
                "--r",
                label="Curvature",
            )
            ax1.plot(
                sig_ela_best_m[i_duct_c_best : i_duct_m_best * 2],
                depth_km[i_duct_c_best : i_duct_m_best * 2],
                "--r",
            )
        else:
            ax1.plot(sig_ela_best, depth_km, "--r", label="Curvature")
        ax1.axvline(sig_y, label="Bounding stress")
        ax1.axvline(-sig_y)
        ax1.set_xlim(
            1.5 * np.min([d_sig_tmp1_best.min(), d_sig_tab2_best.min()]),
            1.5 * np.max([d_sig_tmp0_best.max(), d_sig_tab1_best.max()]),
        )
        ax1.set_xlabel("Sress (Pa)")
        ax1.legend(loc="lower left")
        ax1.set_title(
            "Elastic thickness %i (km), mechanical thickness %i (km)\n"
            % (Te / 1e3, Tm / 1e3)
            + "crust, mantle, and surface heat flows %.2f, %.2f, %.2f (mW m-2)"
            % (F_c_best, F_m_best, F_s_best)
        )

        ax2.plot(T_z_best, depth_km)
        ax2.set_xlabel("Temperature (K)")
        for ax in [ax1, ax2]:
            ax.invert_yaxis()
            ax.set_ylabel("Depth (km)")
        plt.show()
    else:
        ax1 = None
        ax2 = None

    return (
        F_c_best,
        F_m_best,
        F_s_best,
        F_c_best / k_crust,
        F_m_best / k_mantle,
        F_s_best / k_avg,
        Tm,
        d_sig_tmp1_best,
        d_sig_tmp0_best,
        ax1,
        ax2,
    )


def Conversion_Tprofile_Te(
    T_profile,
    z_profile,
    Q_crust,
    A_crust,
    n_crust,
    Q_mantle,
    A_mantle,
    n_mantle,
    g_surf,
    rhobar,
    R_mpr,
    Tc,
    rhoc,
    rhom,
    eps,
    young,
    poisson,
    sig_y,
    R_gas,
    K_curv=1e-7,
    Te_min=5,
    Te_max=300,
    Te_step=1,
    quiet=True,
    plot=False,
    figsize=mpl.rcParams["figure.figsize"],
):
    """
    Determine the elastic thicknesses from the input
    temperature profile and rheology.

    Returns
    -------
    Te : float
       Elastic thickness of the lithosphere (m).
    d_sig_tmp1_best : array, size(max_depth / step_depth)
       The best-fitting YSE compression (Pa).
    d_sig_tmp0_best : array, size(max_depth / step_depth)
       The best-fitting YSE tension (Pa).
    ax1 : matplotlib axes object
        If plot is True – Axis where the YSE is plotted.
    ax2 : matplotlib axes object
        If plot is True – Axis where the temperature profile is plotted.

    Parameters
    ----------
    T_profile : array, float
        Temperature profile (K).
    z_profile : array, float
        Depth (downward positive) associated with each temperature (m).
        The profile should contain with equidistant values and start from 0.
    Q_crust : float
       Activation energy for the crust (J mol-1).
    A_crust : float
       Empirical constant for the crust (Pa-n s-1).
    n_crust : float
       Rheological constant for the crust.
    Q_mantle : float
       Activation energy for the mantle (J mol-1).
    A_mantle : float
       Empirical constant for the mantle (Pa-n s-1).
    n_mantle : float
       Rheological constant for the mantle.
    g_surf : float
       Gravitational attraction of the surface (m s-2).
    rhobar : float
       Mean density of the body (kg m-3).
    R_mpr : float
       Mean radius of the body (m).
    Tc : float
       Crustal thickness (m).
    rhoc : float
       Crustal density (kg m-3).
    rhom : float
       Mantle density (kg m-3).
    eps : float
       Strain rate (s-1).
    young : float
       Young's modulus (Pa).
    poisson : float
       Poisson's ratio.
    sig_y : float
       Bounding stress (Pa).
    R_gas : float
       Gas constant (J mol-1 K-1).
    K_curv : float, optional, default = 1e-7
        Plate curvature (m-1).
    Te_min : float, optional, default = 5
       Minimum elastic thickness investigated (km).
    Te_max : float, optional, default = max(300, z_profile)
       Maximum elastic thickness investigated (km).
    Te_step : float, optional, default = 1
       Elastic thickness step to test (km).
    quiet : boolean, optional, default = True
       if True, print various outputs.
    plot : boolean, optional, default = False
       if True, plot the best-fit yield strength envelope
       and temperature profile.
    figsize : tuple, optional, default = mpl.rcParams['figure.figsize']
       Size of the output figure.
    """

    if z_profile[0] != 0:
        raise ValueError(
            "Input depth profile should start from the surface (depth = 0), here starts at %s"
            % (z_profile[0])
        )

    if z_profile.any() < 0:
        raise ValueError("Input depth should all be positive")

    if np.shape(T_profile) != np.shape(z_profile):
        raise ValueError(
            "T_profile and z_profile should have the same shape. Input was %s and %s"
            % (np.shape(T_profile), np.shape(z_profile))
        )

    max_depth = np.max(z_profile) / 1e3

    if not quiet:
        print("Crustal thickness (km): %i" % (Tc / 1e3))
        print("mantle and crustal densities (kg m-3): %i, %i" % (rhom, rhoc))
        print("Bounding stress (MPa): %i" % (sig_y / 1e6))
        print("Strain rate (s-1): %s" % (eps))
        print("Minimum Te investigated is %.2f" % (Te_min))
        print("Maximum Te investigated is %.2f" % (np.min([Te_max, max_depth])))
        print("Te step investigated is %.2f" % (Te_step))
        print("Curvature (m-1): %s" % (K_curv))

    step_depth = np.abs(z_profile[1] - z_profile[0])
    # step of the integration
    d_sig_tab = np.zeros((3, np.size(z_profile)))

    misfit_temp = 1e50
    decoupling_potential = False
    decoupling = False
    prints_weak = False
    prints_decoup = False
    duct_fail_c = False
    duct_fail_m = False
    i_crust_thick = np.argmin((z_profile - Tc) ** 2)

    # Brittle Strength
    d_sig_tab[0, :], d_sig_tab[1, :] = Brittle_Strength(
        Tc, R_mpr, g_surf, rhom, rhoc, rhobar, z_profile
    )

    # Ductile Strength & Decoupling
    for z, T, iter in zip(z_profile, T_profile, range(len(T_profile))):
        # Ductile
        # Temperature
        if z <= Tc:
            # Onset of ductile deformation
            ductile_sig = (eps / A_crust) ** (1.0 / float(n_crust)) * np.exp(
                Q_crust / (n_crust * R_gas * T)
            )
            d_sig_tab[0, iter] = min(ductile_sig, d_sig_tab[0, iter])
            d_sig_tab[1, iter] = -min(ductile_sig, abs(d_sig_tab[1, iter]))

            # Crustal ductile failure
            if (
                d_sig_tab[0, iter] == ductile_sig or d_sig_tab[1, iter] == -ductile_sig
            ) and not duct_fail_c:
                duct_fail_c = True
                i_duct_c = iter
        else:

            # Onset of ductile deformation
            ductile_sig = (eps / A_mantle) ** (1.0 / float(n_mantle)) * np.exp(
                Q_mantle / (n_mantle * R_gas * T)
            )
            d_sig_tab[0, iter] = min(ductile_sig, d_sig_tab[0, iter])
            d_sig_tab[1, iter] = -min(ductile_sig, abs(d_sig_tab[1, iter]))

            # Mantle ductile failure
            if (
                d_sig_tab[0, iter] == ductile_sig or d_sig_tab[1, iter] == -ductile_sig
            ) and not duct_fail_m:
                duct_fail_m = True
                i_duct_m = iter

    # Bounding stress limit
    if duct_fail_m:
        depth_m0 = np.where(abs(d_sig_tab[0, i_duct_m:]) <= sig_y)
        depth_m1 = np.where(abs(d_sig_tab[1, i_duct_m:]) <= sig_y)
        d_sig_tab[0, i_duct_m:][depth_m0] = 0.0
        d_sig_tab[1, i_duct_m:][depth_m1] = 0.0
    if duct_fail_c:
        depth_c0 = np.where(abs(d_sig_tab[0, i_duct_c:i_crust_thick]) <= sig_y)
        depth_c1 = np.where(abs(d_sig_tab[1, i_duct_c:i_crust_thick]) <= sig_y)
        d_sig_tab[0, i_duct_c:i_crust_thick][depth_c0] = 0.0
        d_sig_tab[1, i_duct_c:i_crust_thick][depth_c1] = 0.0

    M_real, d_sig_tab1, d_sig_tab2, sig_ela, _ = Curv_Moment(
        d_sig_tab[0, :],
        d_sig_tab[1, :],
        K_curv,
        young,
        poisson,
        z_profile.copy(),
        step_depth,
        sig_y,
    )

    if duct_fail_c and not quiet and not prints_weak:
        print("## Apparition of a weak ductile upper layer ##")
        prints_weak = True

    ### Determine whether the rheology is decoupled.
    if (
        duct_fail_c or duct_fail_m and not decoupling_potential
    ):  # Ductile failure required for decoupling
        i_duct = i_duct_c if duct_fail_c else i_duct_m
        if (
            d_sig_tab[0, i_crust_thick:] > 0
        ).any() and (  # Do we have strength below the crust
            (
                sig_ela[i_duct:][depth_m0] > d_sig_tab[0, i_duct:][depth_m0]
            ).any()  # Is the compressional elastic strength larger than ductile failure
            or (
                sig_ela[i_duct:][depth_m1] < d_sig_tab[1, i_duct:][depth_m1]
            ).any()  # Is the extensional elastic strength larger than ductile failure
        ):
            decoupling_potential = True

    if decoupling_potential:
        args_M_real = (K_curv, young, poisson, z_profile.copy(), step_depth, sig_y)

        # Separate all decoupled layers
        decoup = False
        decoup_first = True
        M_real = 0.0
        d_sig_tab1_decoup = np.zeros_like(d_sig_tab[0, :])
        d_sig_tab2_decoup = np.zeros_like(d_sig_tab[0, :])
        i_decoup = []
        i_recoup = []
        for iter in range(len(z_profile)):
            if iter < i_duct:
                continue
            if (
                abs(d_sig_tab[0, iter]) < sig_y or abs(d_sig_tab[1, iter]) < sig_y
            ) and not decoup:
                d_sig_tab_tmp0 = d_sig_tab[0, :].copy()
                d_sig_tab_tmp1 = d_sig_tab[1, :].copy()
                if decoup_first:
                    d_sig_tab_tmp0[iter:] = 0.0
                    d_sig_tab_tmp1[iter:] = 0.0
                else:
                    d_sig_tab_tmp0[: i_recoup[-1]] = 0.0
                    d_sig_tab_tmp1[: i_recoup[-1]] = 0.0
                    d_sig_tab_tmp0[iter:] = 0.0
                    d_sig_tab_tmp1[iter:] = 0.0
                decoup = True
                (
                    M_real_tmp,
                    d_sig_tab1_tmp,
                    d_sig_tab2_tmp,
                    sig_ela_tmp,
                    _,
                ) = Curv_Moment(d_sig_tab_tmp0, d_sig_tab_tmp1, *args_M_real)
                M_real += M_real_tmp
                d_sig_tab1_decoup += d_sig_tab1_tmp
                d_sig_tab2_decoup += d_sig_tab2_tmp
                if decoup_first:
                    sig_ela_decoup = np.array(sig_ela_tmp)
                    decoup_first = False
                else:
                    sig_ela_decoup = np.column_stack((sig_ela_decoup, sig_ela_tmp))
                i_decoup.append(iter)
                if not prints_decoup:
                    decoupling = True
                    prints_decoup = True
                    if not quiet:
                        print(
                            "## Beginning of decoupling between the crust and mantle ##"
                        )
            elif (
                decoup
                and abs(d_sig_tab[0, iter]) > sig_y
                and abs(d_sig_tab[1, iter]) > sig_y
            ):
                decoup = False
                i_recoup.append(iter)

    for Te in (
        np.arange(Te_min, np.min([Te_max, max_depth]), step=Te_step) * 1e3
    ):  # Convert to m
        # Integration of the elastic strength of the lithosphere
        d_sig_tab[2, :] = -(K_curv * young) / float(1.0 - poisson**2)  # Pa
        a = -1
        for j in z_profile:
            a = a + 1
            d_sig_tab[2, a] = (
                -d_sig_tab[2, a] * (j - float(Te) / 2.0) ** 2
            )  # To integrate the elastic strength
            if j > Te:
                bound4 = a
                break
        M_el = np.trapz(
            d_sig_tab[2, :bound4], dx=step_depth
        )  # The analytic formula for M_el is not used here.
        # as we use trapz for M_real, we use it to compute M_el
        # but in practice M_el integrated and the analytic formula are close
        M_el_analyt = K_curv * young * Te**3 / (12.0 * float(1.0 - poisson**2))
        if not quiet:
            print(
                "Numerical integration of M_el vs exact integration M_el",
                "%.4e" % (M_el),
                "%.4e" % (M_el_analyt),
            )

        if np.abs(M_el - M_el_analyt) / M_el_analyt > 5e2:
            raise ValueError(
                "The difference between the integrated elastic moment"
                + " and analytical is too large (%s)"
                % (100 * np.abs(M_el - M_el_analyt) / M_el_analyt)
                + "Increase the depth step sizes, which is here %i km"
                % (step_depth / 1e3)
            )

        misfit = abs(1.0 - abs(M_real / M_el))
        if misfit <= misfit_temp:
            if not quiet:
                print(
                    " YSE misfit %.2f, M_real/M_el %.2f, elastic thickness (km) %.2f "
                    % (misfit, abs(M_real / M_el), Te / 1e3)
                )
            misfit_temp = misfit
            Te_best = Te - 1e3
            d_sig_tmp0_best = d_sig_tab[0, :].copy()
            d_sig_tmp1_best = d_sig_tab[1, :].copy()
            if decoupling:
                sig_ela_best = sig_ela_decoup
                i_decoup_best = i_decoup
                i_recoup_best = i_recoup
                d_sig_tab1_best = d_sig_tab1_decoup
                d_sig_tab2_best = d_sig_tab2_decoup
                decoupling_better = True
            else:
                d_sig_tab1_best = d_sig_tab1
                d_sig_tab2_best = d_sig_tab2
                sig_ela_best = sig_ela
                decoupling_better = False
            find_Tm = np.where(
                np.logical_and(abs(d_sig_tmp0_best) <= sig_y, z_profile > Te_best)
            )
            if np.size(find_Tm) == 0:
                raise ValueError(
                    " Mechanical thickness larger than maximum input depth\n"
                    + "User should increase the max_depth parameter, value is: %i km"
                    % (max_depth)
                )
            else:
                if decoupling_better:
                    Tm1 = 0
                    for ind, i_d in enumerate(i_decoup_best):
                        if ind == 0:
                            Tm2 = z_profile[i_d] ** 3
                        else:
                            Tm2 = (
                                z_profile[i_d] - z_profile[i_recoup_best[ind - 1]]
                            ) ** 3
                        Tm1 += Tm2  # Burov & Diament 1995
                    Tm = Tm1 ** (1 / 3)
                else:
                    Tm = z_profile[find_Tm][0]

        else:
            if not quiet:
                print(
                    " YSE misfit is worse %.2f, M_real/M_el %.2f, elastic thickness %.2f "
                    % (misfit, abs(M_real / M_el), Te / 1e3)
                )

    if not quiet:
        print("Best-fit elastic thickness (km) %.2f " % (Te_best / 1e3))

    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        depth_km = z_profile / 1e3
        ax1.plot(d_sig_tmp0_best, depth_km, "purple")
        ax1.plot(d_sig_tmp1_best, depth_km, "k")
        ax1.plot(d_sig_tab1_best, depth_km, color="orange", label="Integrated tension")
        ax1.plot(d_sig_tab2_best, depth_km, "b", label="Integrated compression")

        if decoupling_better:
            for sigs, i_d, ind in zip(
                sig_ela_best.T if len(i_decoup_best) > 1 else [sig_ela_best.T],
                i_decoup_best,
                range(len(i_decoup_best)),
            ):
                if ind == 0:
                    ax1.plot(
                        sigs[:i_d],
                        depth_km[:i_d],
                        "--r",
                        label="Curvature",
                    )
                else:
                    ax1.plot(
                        sigs[i_recoup_best[ind - 1] : i_d],
                        depth_km[i_recoup_best[ind - 1] : i_d],
                        "--r",
                    )
        else:
            ax1.plot(sig_ela_best, depth_km, "--r", label="Curvature")
        ax1.axvline(sig_y, label="Bounding stress")
        ax1.axvline(-sig_y)
        ax1.set_xlim(
            1.5 * np.min([d_sig_tmp1_best.min(), d_sig_tab2_best.min()]),
            1.5 * np.max([d_sig_tmp0_best.max(), d_sig_tab1_best.max()]),
        )
        ax1.set_xlabel("Sress (Pa)")
        ax1.legend(loc="lower left")
        ax1.set_title(
            "Best-fit elastic thickness %i km, mechanical thickness %i km"
            % (Te_best / 1e3, Tm / 1e3)
        )

        ax2.plot(T_profile, z_profile / 1e3)
        ax2.set_xlabel("Temperature (K)")
        for ax in [ax1, ax2]:
            ax.invert_yaxis()
            ax.set_ylabel("Depth (km)")
        ax2.set_title("Input temperature profile")
        plt.show()
    else:
        ax1 = None
        ax2 = None

    return (Te_best, d_sig_tmp1_best, d_sig_tmp0_best, ax1, ax2)


def Brittle_Strength(Tc, R_mpr, g_surf, rhom, rhoc, rhobar, z_profile):

    """
    Compute the Brittle strength of the lithosphere based on the
    approach of Mueller & Phillips (1995).

    Returns
    -------
    d_sig_tab1 : array size(depth)
       Brittle yield strength envelope in compression (Pa).
    d_sig_tab2 : array size(depth)
       Brittle yield strength envelope in tension (Pa).

    Parameters
    ----------
    Tc : float
       Crustal thickness (m).
    R_mpr : float
       Mean radius of the body (m).
    g_surf : float
       Gravitational attraction of the surface (m s-2).
    rhoc : float
       Crustal density (kg m-3).
    rhom : float
       Mantle density (kg m-3).
    rhobar : float
       Mean density of the body (kg m-3).
    z_profile : array, float
       Depth array (downward positive).
    """

    # gravitational attraction at the moho at the moho
    g_crust = (
        g_surf
        * (1.0 + (((R_mpr - Tc) / R_mpr) ** 3 - 1) * rhoc / rhobar)
        / ((R_mpr - Tc) / R_mpr) ** 2
    )

    z_crust_ix = z_profile <= Tc
    z_mantle_idx = z_profile > Tc
    z_crust = z_profile[z_crust_ix]
    z_mantle = z_profile[z_mantle_idx]

    g_depth = np.hstack(
        (
            g_surf
            * (1.0 + (((R_mpr - z_crust) / R_mpr) ** 3 - 1) * rhoc / rhobar)
            / ((R_mpr - z_crust) / R_mpr) ** 2,  # crust
            g_surf
            * (1.0 + (((R_mpr - z_mantle) / R_mpr) ** 3 - 1) * rhom / rhobar)
            / ((R_mpr - z_mantle) / R_mpr) ** 2,
        )
    )  # mantle

    sig_v = (
        np.hstack(
            (
                rhoc * g_depth[z_crust_ix] * z_crust,  # crust
                rhoc * g_crust * Tc + rhom * g_depth[z_mantle_idx] * (z_mantle - Tc),
            )  # mantle
        )
        * 1e-6
    )

    Tension = (
        np.hstack((0.786 * sig_v[sig_v <= 529.9], 56.7 + 0.679 * sig_v[sig_v > 529.9]))
        * 1e6
    )
    Compression = (
        np.hstack((-3.68 * sig_v[sig_v <= 113.2], -176.6 - 2.12 * sig_v[sig_v > 113.2]))
        * 1e6
    )

    return Tension, Compression
