"""
Convert the elastic thickness into
a heat flow given the input rheology
"""

import numpy as np
import matplotlib.pyplot as plt


def Curv_Moment(
    HF,
    d_sig_tmp0,
    d_sig_tmp1,
    K_curv,
    young,
    poisson,
    depth,
    step_depth,
    sig_y,
    plot_YSE=False,
):
    """
    Determine the bending moment given the input yield
    strength envelope and curvature

    Returns
    -------
    Mx : float
       The bending moment of the plate.
    d_sig_tab1 : array size(depth)
       Integrated part of the yield strength envelope in compression (Pa).
    d_sig_tab2 : array size(depth)
       Integrated part of the yield strength envelope in tension (Pa).
    sigma_ela : array size(depth)
       Integrated elastic part of the yield strength envelope (Pa).

    Parameters
    ----------
    HF : float
       The input mantle heat flow (W m-2).
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
    plot_YSE : boolean, optional, default = False
       If true, plot the yield strength envelope.
    """
    NetAxial = 0
    neutralfib = 0
    a = -1
    first = 0.0
    for i in depth:
        a += 1
        # Find neutral fiber
        sigma_ela = (-young * K_curv * (depth - i)) / (1.0 - poisson ** 2)
        d_sig_tab1 = np.minimum(sigma_ela[:], d_sig_tmp0[:])
        d_sig_tab2 = np.maximum(sigma_ela[:], d_sig_tmp1[:])
        if K_curv >= 0:
            # elastic core limit & neutral fiber
            bound2 = np.min(np.where(sigma_ela == 0))
            d_sig_tab1[bound2:] = 0.0
            d_sig_tab2[:bound2] = 0.0

            NetAxial = np.trapz(d_sig_tab1[:bound2], dx=step_depth) + np.trapz(
                d_sig_tab2[bound2:], dx=step_depth
            )
        else:
            # elastic core limit & neutral fiber
            bound2 = np.min(np.where(sigma_ela == 0))
            d_sig_tab1[:bound2] = 0.0
            d_sig_tab2[bound2:] = 0.0

            NetAxial = np.trapz(d_sig_tab2[:bound2], dx=step_depth) + np.trapz(
                d_sig_tab1[bound2:], dx=step_depth
            )

        if first == 0 and NetAxial != 0:
            tempAx = NetAxial
            first = 1
        if first == 1 and tempAx < 0 and NetAxial > 0:
            neutralfib = i
            break
        elif first == 1 and tempAx > 0 and NetAxial < 0:
            neutralfib = i
            break

    d_sig_tab3 = np.zeros(np.shape(d_sig_tab1))
    d_sig_tab4 = np.zeros(np.shape(d_sig_tab1))

    a = -1  # Update YSE with elastic limit
    for z in depth:
        a += 1
        d_sig_tab3[a] = d_sig_tab1[a] * (z - neutralfib / 2.0)
        d_sig_tab4[a] = d_sig_tab2[a] * (z - neutralfib / 2.0)

    Mx = np.trapz(d_sig_tab3, dx=step_depth) + np.trapz(d_sig_tab4, dx=step_depth)

    if plot_YSE:  # Show YSE
        depth_km = depth / 1e3
        plt.plot(d_sig_tmp0, depth_km, "purple")
        plt.plot(d_sig_tmp1, depth_km, "k")
        plt.plot(d_sig_tab1, depth_km, color="orange", label="Integrated tension")
        plt.plot(d_sig_tab2, depth_km, "b", label="Integrated compression")
        plt.plot(sigma_ela, depth_km, "--r", label="Curvature")
        plt.axvline(sig_y, label="Bounding stress")
        plt.axvline(-sig_y)
        plt.xlim(
            1.5 * np.min([d_sig_tmp1.min(), d_sig_tab2.min()]),
            1.5 * np.max([d_sig_tmp0.max(), d_sig_tmp1.max()]),
        )
        plt.ylabel("Depth (km)")
        plt.xlabel("Sress (Pa)")
        plt.legend(loc="lower left")
        plt.gca().invert_yaxis()
        plt.title("Mantle heat flow %s (mW m-2)" % (HF * 1e3))
        plt.show()

    return Mx, d_sig_tab1, d_sig_tab2, sigma_ela


def Conversion_Te_HF(
    Q_crust,
    A_crust,
    n_crust,
    Q_mantle,
    A_mantle,
    n_mantle,
    g,
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
    Tm : float
       Mechanical thickness of the lithosphere (m).

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
    g : float
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
    step_depth : integer, optional, default = Te / 1000
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
       tested heat flows."""

    if quiet is False:
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
        Tmax = (Q_mantle / R_gas) / (np.log(A_mantle * (sig_y ** n_mantle) / eps))
        F = 1e3 * k_avg * (Tmax - Ts) / Te
    else:  # Crustal rheology
        Tmax = (Q_crust / R_gas) / (np.log(A_crust * (sig_y ** n_crust) / eps))
        F = 1e3 * k_crust * (Tmax - Ts) / Te

    if quiet is False:
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
    depths_arr = np.array(range(1, max_depth, step_depth))

    d_sig_tab = np.zeros((3, np.size(depths_arr)))
    T_z_tab = np.zeros((np.size(depths_arr)))

    # Integration of the elastic strength of the lithosphere
    d_sig_tab[2, :] = -(K_curv * young) / float(1.0 - poisson ** 2)  # Pa
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
    M_el_analyt = K_curv * young * Te ** 3 / (12.0 * float(1.0 - poisson ** 2))
    if quiet is False:
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

    g_crust = (
        g
        * (1.0 + (((R_mpr - Tc) / R_mpr) ** 3 - 1) * rhoc / rhobar)
        / ((R_mpr - Tc) / R_mpr) ** 2
    )  # gravitational attraction at the moho at the moho

    misfit_temp = 1e50
    prints = False
    for HF in HF_arr:  # HF : mw/m2
        duct_fail_c = False
        duct_fail_m = False
        F = float(HF / 1e3)
        d_sig_tab[2, :] = -(K_curv * young) / float(
            1.0 - poisson ** 2
        )  # Strength of the elastic core (without * (z-zn)) Pa
        iter = -1
        for z in depths_arr:
            iter += 1
            if z <= Tc:
                g_depth = (
                    g
                    * (1.0 + (((R_mpr - z) / R_mpr) ** 3 - 1) * rhoc / rhobar)
                    / ((R_mpr - z) / R_mpr) ** 2
                )
            else:
                g_depth = (
                    g
                    * (1.0 + (((R_mpr - z) / R_mpr) ** 3 - 1) * rhom / rhobar)
                    / ((R_mpr - z) / R_mpr) ** 2
                )

            # lithostatic pressure
            if z <= Tc:
                sig_v = rhoc * g_depth * z
            else:
                sig_v = rhoc * g_crust * Tc + rhom * g_depth * (z - Tc)
            sig_v *= 1e-6  # MPa
            # Tension
            if sig_v <= 529.9:  # Mueller & Phillips 1995
                d_sig_tab[0, iter] = float(0.786 * sig_v)
            elif sig_v > 529.9:
                d_sig_tab[0, iter] = float(56.7 + 0.679 * sig_v)
            # Compression
            if sig_v <= 113.2:
                d_sig_tab[1, iter] = float(-3.68 * sig_v)
            elif sig_v > 113.2:
                d_sig_tab[1, iter] = float(-176.6 - 2.12 * sig_v)

            # Ductile
            # Temperature
            if z <= Tc:
                Fs = F + rhoc * H_c * z
                T_z = Ts + Fs * z / k_crust - rhoc * H_c * z ** 2 / (2 * k_crust)

                d_sig_tab[0, iter] = d_sig_tab[0, iter] * 1e6  # Pa
                d_sig_tab[1, iter] = d_sig_tab[1, iter] * 1e6

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
                ) and duct_fail_c is False:
                    duct_fail_c = True
                    iter_duct_c = iter
                    iter_crust_thick = np.argmin((depths_arr - Tc) ** 2)
            else:
                T_z_c = Ts + Fs * Tc / k_crust - rhoc * H_c * Tc ** 2 / (2 * k_crust)
                T_z = T_z_c + F / k_mantle * (z - Tc)

                d_sig_tab[0, iter] = d_sig_tab[0, iter] * 1e6  # Pa
                d_sig_tab[1, iter] = d_sig_tab[1, iter] * 1e6

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
                ) and duct_fail_m is False:
                    duct_fail_m = True
                    iter_duct_m = iter

            T_z_tab[iter] = T_z

        # Bounding stress limit
        if duct_fail_m:
            d_sig_tab[0, iter_duct_m:][
                np.where(abs(d_sig_tab[0, iter_duct_m:]) <= sig_y)
            ] = 0.0
            d_sig_tab[1, iter_duct_m:][
                np.where(abs(d_sig_tab[1, iter_duct_m:]) <= sig_y)
            ] = 0.0
        if duct_fail_c:
            d_sig_tab[0, iter_duct_c:iter_crust_thick][
                np.where(abs(d_sig_tab[0, iter_duct_c:iter_crust_thick]) <= sig_y)
            ] = 0.0
            d_sig_tab[1, iter_duct_c:iter_crust_thick][
                np.where(abs(d_sig_tab[1, iter_duct_c:iter_crust_thick]) <= sig_y)
            ] = 0.0
            if prints is False and d_sig_tab[0, iter_crust_thick + 1] > 0:
                if quiet is False:
                    print("## Beginning of decoupling between the crust and mantle ##")
                prints = True

        M_real, d_sig_tab1, d_sig_tab2, sigma_ela = Curv_Moment(
            F,
            d_sig_tab[0, :].copy(),
            d_sig_tab[1, :].copy(),
            K_curv,
            young,
            poisson,
            depths_arr,
            step_depth,
            sig_y,
            plot_YSE=plot_YSE,
        )

        misfit = abs(1.0 - abs(M_real / M_el))
        if misfit < misfit_temp:
            if quiet is False:
                print(
                    " YSE misfit %.2f, M_real/M_el %.2f, crust, mantle, and surface heat flows (mW m-2) %.2f, %.2f, %.2f "
                    % (misfit, abs(M_real / M_el), (Fs - F) * 1e3, F * 1e3, Fs * 1e3)
                )
            misfit_temp = misfit
            F_m_best = F * 1e3
            F_s_best = Fs * 1e3
            d_sig_tmp0_best = d_sig_tab[0, :].copy()
            d_sig_tmp1_best = d_sig_tab[1, :].copy()
            d_sig_tab1_best = d_sig_tab1
            d_sig_tab2_best = d_sig_tab2
            sigma_ela_best = sigma_ela
            T_z_best = T_z_tab
        else:
            if quiet is False:
                print(
                    " YSE misfit is worse %.2f, M_real/M_el %.2f, crust, mantle, and surface heat flows (mW m-2) %.2f, %.2f, %.2f "
                    % (misfit, abs(M_real / M_el), (Fs - F) * 1e3, F * 1e3, Fs * 1e3)
                )

    # Mechanical thickness of the lithosphere
    Tm = (
        depths_arr[
            np.where(np.logical_and(abs(d_sig_tmp0_best) <= sig_y, depths_arr > Te))
        ][0]
        - 1000
    )
    F_c_best = F_s_best - F_m_best

    if quiet is False:
        print(
            "Best-fit crust, mantle, and surface heat flows (mW m-2) %.2f, %.2f, %.2f"
            % (F_c_best, F_m_best, F_s_best)
        )
        print("Mechanical thickness of the lithosphere (km) %i" % (Tm / 1e3))

    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2)
        depth_km = depths_arr / 1e3
        ax1.plot(d_sig_tmp0_best, depth_km, "purple")
        ax1.plot(d_sig_tmp1_best, depth_km, "k")
        ax1.plot(d_sig_tab1_best, depth_km, color="orange", label="Integrated tension")
        ax1.plot(d_sig_tab2_best, depth_km, "b", label="Integrated compression")
        ax1.plot(sigma_ela_best, depth_km, "--r", label="Curvature")
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

    return (
        F_c_best,
        F_m_best,
        F_s_best,
        F_c_best / k_crust,
        F_m_best / k_mantle,
        F_s_best / k_avg,
        Tm,
    )
