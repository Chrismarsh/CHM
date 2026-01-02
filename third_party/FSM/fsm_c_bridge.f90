!//
!// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
!// modular unstructured mesh based approach for hydrological modelling
!// Copyright (C) 2018 Christopher Marsh
!//
!// This file is part of Canadian Hydrological Model.
!//
!// Canadian Hydrological Model is free software: you can redistribute it and/or
!// modify
!// it under the terms of the GNU General Public License as published by
!// the Free Software Foundation, either version 3 of the License, or
!// (at your option) any later version.
!//
!// Canadian Hydrological Model is distributed in the hope that it will be useful,
!// but WITHOUT ANY WARRANTY; without even the implied warranty of
!// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!// GNU General Public License for more details.
!//
!// You should have received a copy of the GNU General Public License
!// along with Canadian Hydrological Model.  If not, see
!// <http://www.gnu.org/licenses/>.
!//

module fsm_c_bridge
  use iso_c_binding
  implicit none

  interface
     subroutine allocate()
     end subroutine allocate

     subroutine FSM2_TIMESTEP(dt, elev, zT, zU, &
                              LW, Ps, Qa, Rf, Sdif, Sdir, Sf, Ta, trans, Ua, &
                              alb0, hveg, VAI, rhod, &
                              albs, Tsrf, Dsnw, Nsnow, Qcan, Rgrn, Sice, Sliq, &
                              Sveg, Tcan, Tsnow, Tsoil, Tveg, Vsmc, &
                              H, LE, LWout, LWsub, Melt, Roff, snd, snw, subl, svg, &
                              SWout, SWsub, Usub, Wflx)
       real, intent(inout) :: dt
       real, intent(in) :: elev
       real, intent(in) :: zT
       real, intent(in) :: zU
       real, intent(in) :: LW
       real, intent(in) :: Ps
       real, intent(in) :: Qa
       real, intent(inout) :: Rf
       real, intent(in) :: Sdif
       real, intent(in) :: Sdir
       real, intent(inout) :: Sf
       real, intent(in) :: Ta
       real, intent(inout) :: trans
       real, intent(in) :: Ua
       real, intent(in) :: alb0
       real, intent(in) :: hveg
       real, intent(in) :: VAI
       real, intent(in) :: rhod
       real, intent(inout) :: albs
       real, intent(inout) :: Tsrf
       real, intent(inout) :: Dsnw(*)
       integer, intent(inout) :: Nsnow
       real, intent(inout) :: Qcan(*)
       real, intent(inout) :: Rgrn(*)
       real, intent(inout) :: Sice(*)
       real, intent(inout) :: Sliq(*)
       real, intent(inout) :: Sveg(*)
       real, intent(inout) :: Tcan(*)
       real, intent(inout) :: Tsnow(*)
       real, intent(inout) :: Tsoil(*)
       real, intent(inout) :: Tveg(*)
       real, intent(inout) :: Vsmc(*)
       real, intent(out) :: H
       real, intent(out) :: LE
       real, intent(out) :: LWout
       real, intent(out) :: LWsub
       real, intent(out) :: Melt
       real, intent(out) :: Roff
       real, intent(out) :: snd
       real, intent(out) :: snw
       real, intent(out) :: subl
       real, intent(out) :: svg
       real, intent(out) :: SWout
       real, intent(out) :: SWsub
       real, intent(out) :: Usub
       real, intent(out) :: Wflx(*)
     end subroutine FSM2_TIMESTEP
  end interface

contains

  subroutine fsm_allocate() bind(C, name="fsm_allocate")
    call allocate()
  end subroutine fsm_allocate

  subroutine fsm_layers_config(fvg1_in, zsub_in, ncnpy_in, nsmax_in, nsoil_in) bind(C, name="fsm_layers_config")
    use layers, only: fvg1, zsub, Ncnpy, Nsmax, Nsoil
    real(c_float), value :: fvg1_in
    real(c_float), value :: zsub_in
    integer(c_int), value :: ncnpy_in
    integer(c_int), value :: nsmax_in
    integer(c_int), value :: nsoil_in

    fvg1 = fvg1_in
    zsub = zsub_in
    Ncnpy = ncnpy_in
    Nsmax = nsmax_in
    Nsoil = nsoil_in
  end subroutine fsm_layers_config

  subroutine fsm_soilprops_config(b_in, hcap_in, hcon_in, sathh_in, vcrit_in, vsat_in) bind(C, name="fsm_soilprops_config")
    use soilprops, only: b, hcap_soil, hcon_soil, sathh, Vcrit, Vsat
    real(c_float), value :: b_in
    real(c_float), value :: hcap_in
    real(c_float), value :: hcon_in
    real(c_float), value :: sathh_in
    real(c_float), value :: vcrit_in
    real(c_float), value :: vsat_in

    b = b_in
    hcap_soil = hcap_in
    hcon_soil = hcon_in
    sathh = sathh_in
    Vcrit = vcrit_in
    Vsat = vsat_in
  end subroutine fsm_soilprops_config

  function fsm_constants_e0() result(val) bind(C, name="fsm_constants_e0")
    use constants, only: e0
    real(c_float) :: val

    val = e0
  end function fsm_constants_e0

  function fsm_constants_eps() result(val) bind(C, name="fsm_constants_eps")
    use constants, only: eps
    real(c_float) :: val

    val = eps
  end function fsm_constants_eps

  function fsm_parameters_rgr0() result(val) bind(C, name="fsm_parameters_rgr0")
    use parameters, only: rgr0
    real(c_float) :: val

    val = rgr0
  end function fsm_parameters_rgr0

  subroutine fsm2_timestep_c(dt, elev, zT, zU, &
                           LW, Ps, Qa, Rf, Sdif, Sdir, Sf, Ta, trans, Ua, &
                           alb0, hveg, VAI, rhod, &
                           albs, Tsrf, Dsnw, Nsnow, Qcan, Rgrn, Sice, Sliq, &
                           Sveg, Tcan, Tsnow, Tsoil, Tveg, Vsmc, &
                           H, LE, LWout, LWsub, Melt, Roff, snd, snw, subl, svg, &
                           SWout, SWsub, Usub, Wflx) bind(C, name="fsm2_timestep")
    real(c_float), intent(inout) :: dt
    real(c_float), intent(in) :: elev
    real(c_float), intent(in) :: zT
    real(c_float), intent(in) :: zU
    real(c_float), intent(in) :: LW
    real(c_float), intent(in) :: Ps
    real(c_float), intent(in) :: Qa
    real(c_float), intent(inout) :: Rf
    real(c_float), intent(in) :: Sdif
    real(c_float), intent(in) :: Sdir
    real(c_float), intent(inout) :: Sf
    real(c_float), intent(in) :: Ta
    real(c_float), intent(inout) :: trans
    real(c_float), intent(in) :: Ua
    real(c_float), intent(in) :: alb0
    real(c_float), intent(in) :: hveg
    real(c_float), intent(in) :: VAI
    real(c_float), intent(in) :: rhod
    real(c_float), intent(inout) :: albs
    real(c_float), intent(inout) :: Tsrf
    real(c_float), intent(inout) :: Dsnw(*)
    integer(c_int), intent(inout) :: Nsnow
    real(c_float), intent(inout) :: Qcan(*)
    real(c_float), intent(inout) :: Rgrn(*)
    real(c_float), intent(inout) :: Sice(*)
    real(c_float), intent(inout) :: Sliq(*)
    real(c_float), intent(inout) :: Sveg(*)
    real(c_float), intent(inout) :: Tcan(*)
    real(c_float), intent(inout) :: Tsnow(*)
    real(c_float), intent(inout) :: Tsoil(*)
    real(c_float), intent(inout) :: Tveg(*)
    real(c_float), intent(inout) :: Vsmc(*)
    real(c_float), intent(out) :: H
    real(c_float), intent(out) :: LE
    real(c_float), intent(out) :: LWout
    real(c_float), intent(out) :: LWsub
    real(c_float), intent(out) :: Melt
    real(c_float), intent(out) :: Roff
    real(c_float), intent(out) :: snd
    real(c_float), intent(out) :: snw
    real(c_float), intent(out) :: subl
    real(c_float), intent(out) :: svg
    real(c_float), intent(out) :: SWout
    real(c_float), intent(out) :: SWsub
    real(c_float), intent(out) :: Usub
    real(c_float), intent(out) :: Wflx(*)

    call FSM2_TIMESTEP(dt, elev, zT, zU, &
                       LW, Ps, Qa, Rf, Sdif, Sdir, Sf, Ta, trans, Ua, &
                       alb0, hveg, VAI, rhod, &
                       albs, Tsrf, Dsnw, Nsnow, Qcan, Rgrn, Sice, Sliq, &
                       Sveg, Tcan, Tsnow, Tsoil, Tveg, Vsmc, &
                       H, LE, LWout, LWsub, Melt, Roff, snd, snw, subl, svg, &
                       SWout, SWsub, Usub, Wflx)
  end subroutine fsm2_timestep_c

end module fsm_c_bridge
