module chm_crocus_bridge
    use iso_c_binding, only: c_float, c_int
    implicit none
    private
    public :: chm_crocus_initialize, chm_crocus_step

    ! Keep these interoperable types in the same order as crocus_bridge.hpp.
    type, bind(C) :: forcing
        integer(c_int) :: year, month, day, seconds
        real(c_float) :: dt, air_k, pressure_pa, specific_humidity, wind, sw, lw
        real(c_float) :: snowfall, rainfall, cos_slope, solar_zenith
        real(c_float) :: transported_mass, transported_density
    end type
    type, bind(C) :: boundary
        real(c_float) :: ground_k, conductivity, ground_layer_m, ground_albedo
    end type
    type, bind(C) :: fluxes
        real(c_float) :: albedo, sensible, latent, lw_out, runoff, sublimation, surface_melt, basal_heat
    end type
contains
    subroutine chm_crocus_initialize() bind(C)
        use modi_ini_csts
        use modi_ini_surf_csts
        use modd_prep_snow, only: nimpur
        nimpur = 0
        call ini_csts()
        call ini_surf_csts()
    end subroutine

    subroutine chm_crocus_step(n, state, met, ground, pending_snow, out) bind(C)
        use modd_csts, only: xrd, xrv, xcpd, xp00, xg, xstefan, xlstt
        use modd_type_date_surf, only: date_time
        use modd_snow_metamo, only: xuepsi_smp, xuepsi
        use modi_snowcro
        use sfc_options, only: hsnowres
        integer(c_int), value :: n
        real(c_float), intent(inout) :: state(1,n,7), pending_snow
        type(forcing), intent(in) :: met
        type(boundary), intent(in) :: ground
        type(fluxes), intent(inout) :: out
        type(date_time) :: today
        real :: snow_rate, rho_air, exner, exner_air, wind(1), resistance(1), blowing(1,4)
        real :: liquid(1,n), temperature(1,n), thickness(1,n), impurity(1,n,0), deposition(1,0)
        real :: albedo(1), drainage(1), basal(1), evap_correction(1), ground_correction(1)
        real :: short_net(1), short_surface(1), long_net(1), radiation(1), sensible(1), net_heat(1)
        real :: rain_heat(1), latent_ice(1), latent_water(1), evap(1), drift_loss(1), ri(1)
        real :: emissivity(1), drag(1), ustar(1), heat_drag(1), mass_heat(1), qs(1), snow_flux(1), melt(1)
        real :: band_albedo(1,1), band_diffuse(1,1), band_total(1,1)
        real :: incident(1,1)
        real, parameter :: zref = 2.0, roughness = 0.001, thermal_roughness = 0.0001

        ! The upstream fresh-snow routine has a single-precision mass threshold.
        ! Accumulate sub-threshold snowfall instead of losing water or calling with zero layers.
        blowing = 0.
        if (met%transported_mass > 200.0*xuepsi_smp .and. met%transported_mass/met%dt > xuepsi) then
            ! Native Crocus accepts deposited snow separately from atmospheric snowfall.
            blowing(1,:) = [met%transported_mass/met%dt, met%transported_density, 0.0003, 0.9]
        else
            pending_snow = pending_snow + met%transported_mass
        end if
        pending_snow = pending_snow + met%snowfall * met%dt
        snow_rate = 0.
        if (pending_snow > 200.0*xuepsi_smp) then
            snow_rate = pending_snow / met%dt
            pending_snow = 0.
        end if
        if (state(1,1,1) <= 0. .and. snow_rate <= 0. .and. blowing(1,1) <= 0.) then
            out = fluxes(ground%ground_albedo, 0., 0., xstefan*ground%ground_k**4, met%rainfall, 0., 0., 0.)
            return
        end if

        today%tdate%year = met%year
        today%tdate%month = met%month
        today%tdate%day = met%day
        today%time = real(met%seconds)
        rho_air = met%pressure_pa / (xrd*met%air_k*(1.+(xrv/xrd-1.)*met%specific_humidity))
        exner = (met%pressure_pa/xp00)**(xrd/xcpd)
        exner_air = ((met%pressure_pa-rho_air*xg*zref)/xp00)**(xrd/xcpd)
        albedo = out%albedo
        wind = met%wind
        resistance = -1. ! use the native surface-exchange calculation
        incident = 0.
        basal = 0.
        short_net = 0.
        short_surface = 0.
        long_net = 0.
        radiation = 0.
        sensible = 0.
        net_heat = 0.
        rain_heat = 0.
        latent_ice = 0.
        latent_water = 0.
        evap = 0.
        ri = 0.
        drag = 0.
        ustar = 0.
        heat_drag = 0.

        ! The ECCC open-column configuration uses prescribed meteorology (A=0,
        ! B=forcing), and lets Crocus solve the surface energy balance (OMEB=false).
        call snowcro(hsnowres, today, .false., 'NEW', [0.], wind, &
            [0.], [0.], [met%air_k/exner], [met%specific_humidity], &
            state(:,:,1), state(:,:,2), state(:,:,3), albedo, &
            state(:,:,4), state(:,:,5), state(:,:,6), state(:,:,7), impurity, &
            met%dt, [met%pressure_pa], [snow_rate], [0.], [met%rainfall], [1.], resistance, &
            [met%air_k], [ground%ground_k], [met%sw], [met%specific_humidity], wind, wind, [met%lw], [rho_air], &
            [zref], [exner], [exner_air], [met%cos_slope], [0.], &
            [zref], [roughness], [roughness], [thermal_roughness], [ground%ground_albedo], &
            [ground%conductivity], [ground%ground_layer_m], liquid, temperature, thickness, &
            drainage, basal, evap_correction, ground_correction, short_net, short_surface, long_net, &
            radiation, sensible, net_heat, rain_heat, latent_ice, latent_water, evap, drift_loss, ri, &
            emissivity, drag, ustar, heat_drag, mass_heat, qs, [0.], [met%solar_zenith], [0.], [0.], [0.], blowing, &
            'VI13', 'NONE', .false., .false., 'B21', 'B92', .false., incident, incident, &
            band_albedo, band_diffuse, band_total, snow_flux, deposition, deposition, &
            'V12', 'Y81', 'B92', 'B92', 'CST', [0.], .false., .false., .false., .false., .false., &
            [0.], [0.], .false., [60.], .false., melt)

        out%albedo = albedo(1)
        if (state(1,1,1) <= 0.) out%albedo = ground%ground_albedo
        out%sensible = sensible(1)
        out%sublimation = evap(1) - evap_correction(1)
        out%latent = latent_ice(1) + latent_water(1) - xlstt*evap_correction(1)
        out%lw_out = met%lw - (radiation(1) - met%sw*(1.-albedo(1)))
        out%runoff = drainage(1)
        out%surface_melt = melt(1)
        out%basal_heat = basal(1)
    end subroutine
end module
