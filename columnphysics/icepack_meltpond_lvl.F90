!=======================================================================

! Level-ice meltpond parameterization
!
! This meltpond parameterization was developed for use with the delta-
! Eddington radiation scheme, and only affects the radiation budget in
! the model.  That is, although the pond volume is tracked, that liquid
! water is not used elsewhere in the model for mass budgets or other
! physical processes.
!
! authors Elizabeth Hunke (LANL)
!         David Hebert (NRL Stennis)
!         Olivier Lecomte (Univ. Louvain)

      module icepack_meltpond_lvl

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, c10, p01, p5, puny, p333
      use icepack_parameters, only: viscosity_dyn, rhoi, rhos, rhow, Timelt, Tffresh, Lfresh
      use icepack_parameters, only: gravit, depressT, rhofresh, kice, pndaspect, use_smliq_pnd
      use icepack_parameters, only: ktherm, frzpnd, dpscale, hi_min
      use icepack_parameters, only: pndfrbd, pndhyps, pndhead
      use icepack_parameters, only: rhosi, apond_sl
      use icepack_parameters, only: pnd_L_hi0, pnd_L_c02, pnd_L_c0, pnd_L_c1, pnd_L_c2
      use icepack_parameters, only: pnd_a0_hi0, pnd_a0_c03, pnd_a0_c02, pnd_a0_c0, pnd_a0_c1
      use icepack_parameters, only: pnd_y_hi0, pnd_y_c03, pnd_y_c02, pnd_y_c0, pnd_y_c1
      use icepack_parameters, only: pnd_k, pnd_hi_min
      use icepack_tracers,    only: nilyr
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: compute_ponds_lvl, pond_hypsometry, pond_head

!=======================================================================

      contains

!=======================================================================

      subroutine compute_ponds_lvl(dt,                   &
                                   rfrac,  meltt, melts, &
                                   frain,  Tair,  fsurfn,&
                                   dhs,    ffrac,        &
                                   aicen,  vicen, vsnon, &
                                   qicen,  sicen,        &
                                   Tsfcn,  alvl,         &
                                   apnd,   hpnd,  ipnd,  &
                                   meltsliqn, frpndn,    &
                                   rfpndn, ilpndn)

      real (kind=dbl_kind), intent(in) :: &
         dt          ! time step (s)

      real (kind=dbl_kind), intent(in) :: &
         Tsfcn, &    ! surface temperature (C)
         alvl,  &    ! fraction of level ice
         rfrac, &    ! water fraction retained for melt ponds
         meltt, &    ! top melt rate (m/s)
         melts, &    ! snow melt rate (m/s)
         frain, &    ! rainfall rate (kg/m2/s)
         Tair,  &    ! air temperature (K)
         fsurfn,&    ! atm-ice surface heat flux  (W/m2)
         aicen, &    ! ice area fraction
         vicen, &    ! ice volume (m)
         vsnon, &    ! snow volume (m)
         meltsliqn   ! liquid contribution to meltponds in dt (kg/m^2)

      real (kind=dbl_kind), intent(inout) :: &
         apnd, hpnd, ipnd, &
         frpndn, &   ! pond drainage rate due to freeboard constraint (m/step)
         rfpndn, &   ! runoff rate due to rfrac (m/step)
         ilpndn      ! pond loss/gain due to ice lid (m/step)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         qicen, &  ! ice layer enthalpy (J m-3)
         sicen     ! salinity (ppt)

      real (kind=dbl_kind), intent(in) :: &
         dhs       ! depth difference for snow on sea ice and pond ice

      real (kind=dbl_kind), intent(out) :: &
         ffrac     ! fraction of fsurfn over pond used to melt ipond

      ! local temporary variables

      real (kind=dbl_kind) :: &
         volpn, &     ! pond volume per unit area (m)
         dhpondn , &  ! vertical change in pond height (m)
         dvn_temp     ! local variable for change in volume due to rfrac

      real (kind=dbl_kind), dimension (nilyr) :: &
         Tmlt      ! melting temperature (C)

      real (kind=dbl_kind) :: &
         hi                     , & ! ice thickness (m)
         hs                     , & ! snow depth (m)
         dTs                    , & ! surface temperature diff for freeze-up (C)
         Tp                     , & ! pond freezing temperature (C)
         Ts                     , & ! surface air temperature (C)
         apondn                 , & ! local pond area
         hpondn                 , & ! local pond depth (m)
         dvn                    , & ! change in pond volume (m)
         hlid, alid             , & ! refrozen lid thickness, area
         dhlid                  , & ! change in refrozen lid thickness
         bdt                    , & ! 2 kice dT dt / (rhoi Lfresh)
         alvl_tmp               , & ! level ice fraction of ice area
         draft, deltah, pressure_head, perm, drain ! for permeability

      real (kind=dbl_kind), parameter :: &
         Td       = c2          , & ! temperature difference for freeze-up (C)
         rexp     = p01             ! pond contraction scaling
      
      character (len=char_len) :: &
         hypso_type        ! string indicating how to change pond depth-area

      character(len=*),parameter :: subname='(compute_ponds_lvl)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      volpn = hpnd * aicen * alvl * apnd
      ffrac = c0

      !-----------------------------------------------------------------
      ! Identify grid cells where ponds can be
      !-----------------------------------------------------------------

      if (aicen*alvl > puny**2) then

         hi = vicen/aicen
         hs = vsnon/aicen
         alvl_tmp = alvl

         if (hi < hi_min) then

            !--------------------------------------------------------------
            ! Remove ponds on thin ice
            !--------------------------------------------------------------
            apondn = c0
            hpondn = c0
            volpn  = c0
            hlid = c0

         else

            !-----------------------------------------------------------
            ! initialize pond area as fraction of ice
            !-----------------------------------------------------------
            apondn = apnd*alvl_tmp

            !-----------------------------------------------------------
            ! update pond volume
            !-----------------------------------------------------------
            ! add melt water
            if (use_smliq_pnd) then
               dvn = rfrac/rhofresh*(meltt*rhoi &
                   +                 meltsliqn)*aicen
            else
               dvn = rfrac/rhofresh*(meltt*rhoi &
                   +                 melts*rhos &
                   +                 frain*  dt)*aicen
            endif
            ! Track lost meltwater dvn is volume of meltwater (m3/m2) captured
            ! over entire grid cell area. Multiply by (1-rfrac)/rfrac to get
            ! loss over entire area. And divide by aicen to get loss per unit
            ! category area (for consistency with melttn, frpndn, etc)
            rfpndn = dvn * (c1-rfrac) / (rfrac * aicen)
            dvn_temp = dvn

            ! shrink pond volume under freezing conditions
            if (trim(frzpnd) == 'cesm') then
               Tp = Timelt - Td
               dTs = max(Tp - Tsfcn,c0)
               dvn = dvn - volpn * (c1 - exp(rexp*dTs/Tp))

            else
               ! trim(frzpnd) == 'hlid' Stefan approximation
               ! assumes pond is fresh (freezing temperature = 0 C)
               ! and ice grows from existing pond ice
               hlid = ipnd
               if (dvn == c0) then ! freeze pond
                  Ts = Tair - Tffresh
                  if (Ts < c0) then
                     ! if (Ts < -c2) then ! as in meltpond_cesm
                     bdt = -c2*Ts*kice*dt/(rhoi*Lfresh)
                     dhlid = p5*sqrt(bdt)                  ! open water freezing
                     if (hlid > dhlid) dhlid = p5*bdt/hlid ! existing ice
                     dhlid = min(dhlid, hpnd*rhofresh/rhoi)
                     hlid = hlid + dhlid
                  else
                     dhlid = c0 ! to account for surface inversions
                  endif
               else ! convert refrozen pond ice back to water
                  dhlid = max(fsurfn*dt / (rhoi*Lfresh), c0) ! > 0
                  dhlid = -min(dhlid, hlid) ! < 0
                  hlid = max(hlid + dhlid, c0)
                  if (hs - dhs < puny) then ! pond ice is snow-free
                     ffrac = c1 ! fraction of fsurfn over pond used to melt ipond
                     if (fsurfn > puny) &
                          ffrac = min(-dhlid*rhoi*Lfresh/(dt*fsurfn), c1)
                  endif
               endif
               alid = apondn * aicen
               dvn = dvn - dhlid*alid*rhoi/rhofresh
            endif

            ! Track lost/gained meltwater per unit category area from pond 
            ! lid freezing/melting. Note sign flip relative to dvn convention
            ilpndn = (dvn_temp - dvn) / aicen

            !-----------------------------------------------------------
            ! update pond area and depth
            !-----------------------------------------------------------
            dhpondn = c0
            if (trim(pndhyps) == 'none') then
               hypso_type = 'aspect_change'
               call pond_hypsometry(volpn, apondn, hpondn, dvn, alvl_tmp, &
                                    aicen, hypso_type, dhpondn, hi)
            elseif (trim(pndhyps) == 'fixed') then
               hypso_type = 'aspect_fixed'
               call pond_hypsometry(volpn, apondn, hpondn, dvn, alvl_tmp, & 
                                    aicen, hypso_type, dhpondn, hi)
            elseif (trim(pndhyps) == 'sealevel_lin') then
               hypso_type = 'sealevel_lin'
               call pond_hypsometry(volpn, apondn, hpondn, dvn, alvl_tmp, & 
                                    aicen, hypso_type, dhpondn, hi)
            elseif (trim(pndhyps) == 'sealevel_log') then
               hypso_type = 'sealevel_log'
               call pond_hypsometry(volpn, apondn, hpondn, dvn, alvl_tmp, & 
                                    aicen, hypso_type, dhpondn, hi)
            else
               call icepack_warnings_add(subname//" invalid pndhyps option" )
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               if (icepack_warnings_aborted(subname)) return
            endif

            ! limit pond depth to maintain nonnegative freeboard
            if (trim(pndfrbd) == 'floor') then
               dhpondn = min(c0, ((rhow-rhoi)*hi - rhos*hs)/rhofresh - hpondn)
               dvn = dhpondn * aicen * apondn
            elseif (trim(pndfrbd) == 'category') then
               dhpondn = min(c0, ((rhow-rhoi)*hi - rhos*hs)/(rhofresh*apondn) & 
                                 - hpondn)
               dvn = dhpondn * aicen * apondn
            else
               call icepack_warnings_add(subname//" invalid pndfrbd option" )
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               if (icepack_warnings_aborted(subname)) return
            endif
            if (trim(pndhyps) == 'none') then
               hypso_type = 'vertical'
               call pond_hypsometry(volpn, apondn, hpondn, dvn, alvl_tmp, &
                                    aicen, hypso_type, dhpondn, hi)
               frpndn = -dhpondn * apondn
            elseif (trim(pndhyps) == 'fixed') then
               hypso_type = 'aspect_fixed'
               call pond_hypsometry(volpn, apondn, hpondn, dvn, alvl_tmp, & 
                                    aicen, hypso_type, dhpondn, hi)
               ! Meltwater volume change per grid cell area divided by aicen here 
               ! yields the meltwater volume lost averaged over the category area
               ! analogous to how melttn is defined. Note sign flip
               frpndn = - dvn / aicen
            elseif (trim(pndhyps) == 'sealevel') then
               hypso_type = 'sealevel'
               call pond_hypsometry(volpn, apondn, hpondn, dvn, alvl_tmp, & 
                                    aicen, hypso_type, dhpondn, hi)
               frpndn = - dvn / aicen
            elseif (trim(pndhyps) == 'sealevel_log') then
               hypso_type = 'sealevel_log'
               call pond_hypsometry(volpn, apondn, hpondn, dvn, alvl_tmp, & 
                                    aicen, hypso_type, dhpondn, hi)
               frpndn = - dvn / aicen
            else
                  call icepack_warnings_add(subname//" invalid pndhyps option" )
                  call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
                  if (icepack_warnings_aborted(subname)) return
             endif
            
            ! fraction of grid cell covered by ponds
            apondn = apondn * aicen

            volpn = hpondn*apondn
            ! note, this implies that if ponds fully drain or freeze their
            ! depressions cease to exist and the lid ice also ceases to exist
            if (volpn <= c0) then
               volpn = c0
               apondn = c0
               hpondn = c0
               hlid = c0
            endif

            !-----------------------------------------------------------
            ! drainage due to permeability (flushing)
            ! setting dpscale = 0 turns this off
            ! NOTE this uses the initial salinity and melting T profiles
            !-----------------------------------------------------------

            if (ktherm /= 2 .and. hpondn > c0 .and. dpscale > puny) then
               draft = (rhos*hs + rhoi*hi)/rhow + hpondn
               deltah = hpondn + hi - draft
               pressure_head = gravit * rhow * max(deltah, c0)
               Tmlt(:) = -sicen(:) * depressT
               call brine_permeability(qicen, &
                    sicen, Tmlt, perm)
               if (icepack_warnings_aborted(subname)) return
               drain = perm*pressure_head*dt / (viscosity_dyn*hi) * dpscale
               deltah = min(drain, hpondn)
               dvn = -deltah*apondn
               volpn = volpn + dvn
               apondn = max(c0, min(apondn &
                    + 0.5*dvn/(pndaspect*apondn), alvl_tmp*aicen))
               hpondn = c0
               if (apondn > puny) hpondn = volpn/apondn
            endif

         endif

         !-----------------------------------------------------------
         ! Reload tracer array
         !-----------------------------------------------------------

         hpnd = hpondn
         apnd = apondn / (aicen*alvl_tmp)
         if (trim(frzpnd) == 'hlid') ipnd = hlid

      endif

      end subroutine compute_ponds_lvl

!=======================================================================

! determine the liquid fraction of brine in the ice and the permeability

      subroutine brine_permeability(qicen, salin, Tmlt, perm)

      use icepack_therm_shared, only: calculate_Tin_from_qin

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         qicen, &  ! enthalpy for each ice layer (J m-3)
         salin, &  ! salinity (ppt)
         Tmlt      ! melting temperature (C)

      real (kind=dbl_kind), intent(out) :: &
         perm      ! permeability (m^2)

      ! local variables

      real (kind=dbl_kind) ::   &
         Sbr       ! brine salinity

      real (kind=dbl_kind), dimension(nilyr) ::   &
         Tin, &    ! ice temperature (C)
         phi       ! liquid fraction

      integer (kind=int_kind) :: k

      character(len=*),parameter :: subname='(brine_permeability)'

      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Tin(k) = calculate_Tin_from_qin(qicen(k),Tmlt(k))
      enddo

      !-----------------------------------------------------------------
      ! brine salinity and liquid fraction
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Sbr = c1/(1.e-3_dbl_kind - depressT/Tin(k)) ! Notz thesis eq 3.6
         phi(k) = salin(k)/Sbr ! liquid fraction
         if (phi(k) < 0.05) phi(k) = c0 ! impermeable
      enddo

      !-----------------------------------------------------------------
      ! permeability
      !-----------------------------------------------------------------

      perm = 3.0e-8_dbl_kind * (minval(phi))**3

      end subroutine brine_permeability

!=======================================================================

! Manage how pond depth and area changes due to changing water volume

      subroutine pond_hypsometry(volpn, apondn, hpondn,  &
                                 dvn,   alvln,  aicen,   &
                                 type,  dhpondn, hin)

      real (kind=dbl_kind), intent(in) :: &
         hin     , & ! category mean ice thickness
         dvn     , & ! change in pond volume per unit area of grid cell (m)
         alvln   , & ! level ice fraction of category area
         aicen       ! category area fraction
      
      character (len=char_len), intent(in) :: &
         type        ! how to change pond depth and area

      real (kind=dbl_kind), intent(inout) :: &
         dhpondn , & ! change in pond depth, only used if type=='vertical'
         volpn   , & ! pond volume per unit area of the grid cell (m)
         apondn  , & ! pond area fraction of the category (incl. deformed)
         hpondn      ! pond depth (m)
      
      real (kind=dbl_kind) :: &
         temp    , & ! temporary variable for recording changes
         L       , & ! amplitude parameter for logistic hypsometry
         a0      , & ! location parameter for logistic hypsometry
         y       , & ! vertical offset parameter for logistic hypsometry
         Vp      , & ! pond volume per unit area of the category (m)
         V           ! temporary pond volume per unit area of category (m)

      integer (kind=int_kind) :: &
         i           ! counter for Newton's method
      
      real (kind=dbl_kind), parameter :: &
         tol = 0.00001_dbl_kind  ! Tolerance for hypsometry

      integer (kind=int_kind), parameter :: &
         maxi = 100              ! Maximum number of iterations
      
      character(len=*),parameter :: subname='(pond_hypsometry)'

      ! Change pond depth and area
      if (trim(type) == 'vertical') then
         ! Only change depth, not area
         temp = hpondn
         hpondn = max(c0, hpondn + dhpondn)
         volpn = hpondn * aicen * apondn
         if (hpondn == c0) then
            apondn = c0
         endif
         dhpondn = hpondn - temp
      elseif (trim(type) == 'aspect_change') then
         ! update pond volume
         volpn = volpn + dvn
         if (volpn <= c0) then
            volpn = c0
            apondn = c0
            hpondn = c0
         endif
         ! Apply pndaspect ratio to change in depth and area
         if (apondn*aicen > puny) then ! existing ponds
            apondn = max(c0, min(alvln, &
                 apondn + 0.5*dvn/(pndaspect*apondn*aicen)))
            hpondn = c0
            if (apondn > puny) &
                 hpondn = volpn/(apondn*aicen)

         elseif (alvln*aicen > c10*puny) then ! new ponds
            apondn = min (sqrt(volpn/(pndaspect*aicen)), alvln)
            hpondn = pndaspect * apondn

         else           ! melt water runs off deformed ice
            apondn = c0
            hpondn = c0
         endif
         apondn = max(apondn, c0)
      elseif (trim(type) == 'aspect_fixed') then
         ! update pond volume
         volpn = volpn + dvn
         if (volpn <= c0) then
            volpn = c0
            apondn = c0
            hpondn = c0
         endif

         ! compute the pond area and depth
         apondn = min(alvln, sqrt(volpn/(pndaspect*aicen)))
         if (apondn >= alvln) then ! pond fills all available space
            hpondn = volpn / (aicen * apondn)
         else ! pond follows aspect ratio
            hpondn = pndaspect * apondn
         endif
      elseif (trim(type) == 'sealevel_lin') then
         ! update pond volume
         volpn = volpn + dvn
         if (volpn <= c0) then
            volpn = c0
            apondn = c0
            hpondn = c0
         endif
         ! Compute the pond aspect ratio for sea level ponds
         pndaspect = hin*(rhow - rhosi) / &
                     (rhofresh*apond_sl**c2 - c2*rhow*apond_sl + rhow)
         ! compute the pond area and depth
         apondn = min(alvln, sqrt(volpn/(pndaspect*aicen)))
         if (apondn >= alvln) then ! pond fills all available space
            hpondn = volpn / (aicen * apondn)
         else ! pond follows aspect ratio
            hpondn = pndaspect * apondn
         endif
      elseif (trim(type) == 'sealevel_log') then
         ! update pond volume
         volpn = volpn + dvn
         if (volpn <= c0) then
            volpn = c0
            apondn = c0
            hpondn = c0
         else
            ! If the ice is thin enough just use sealevel_lin ponds
            if (hin < pnd_hi_min) then
               ! Compute the pond aspect ratio for sea level ponds
               pndaspect = hin*(rhow - rhosi)/ &
                           (rhofresh*apond_sl**c2 - c2*rhow*apond_sl + rhow)
               ! compute the pond area and depth
               apondn = min(c1, sqrt(volpn/(pndaspect*aicen)))
               if (apondn >= c1) then ! pond fills all available space
                  hpondn = volpn/aicen
               else ! pond follows aspect ratio
                  hpondn = pndaspect*apondn
               endif
            else ! Thick ice, use logistic hypsometry
               ! Calculate hypsometric parameters
               L = pond_log_L(hin)
               a0 = pond_log_a0(hin)
               y = pond_log_y(hin)
               Vp = volpn/aicen
               ! Check if pond will fill entire category
               if (Vp >= (pond_log_e(c1,L,pnd_k,a0,y) &
                         - pond_log_inte(c1,L,pnd_k,a0,y))) then
                  apondn = c1
                  hpondn = Vp
               else
                  ! Use newtons method to solve for apondn
                  apondn = apond_sl ! initial guess for apondn
                  V = pond_log_Vp(apondn,L,pnd_k,a0,y)
                  i = 0
                  do while ((ABS(V - Vp) > tol) .and. (i <= maxi))
                     apondn = apondn - V &
                              /(apondn*pond_log_deda(apondn,L,pnd_k,a0,y))
                     apondn = max(apondn, c0)
                     apondn = min(apondn, c1)
                     V = pond_log_Vp(apondn,L,pnd_k,a0,y)
                     i = i + 1
                  enddo
                  if (apondn <= puny) then
                     apondn = c0
                     hpondn = c0
                  else
                     hpondn = V / apondn
                  endif
               endif
            endif ! hin < pnd_hi_min
         endif ! volpn =< c0
      else
         call icepack_warnings_add(subname//" invalid type option" )
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         if (icepack_warnings_aborted(subname)) return
      endif ! type

      end subroutine pond_hypsometry

!=======================================================================

      function pond_log_L(hi) result(L)

         ! Polynomial approximation of L for logistic pond hypsometry

         real(kind=dbl_kind), intent(in) :: &
            hi    ! ice thickness (m)
         
         real(kind=dbl_kind) :: &
            L     ! amplitude parameter for logistic hypsometry

         character(len=*),parameter :: subname='(pond_log_L)'

         L = pnd_L_c02*(hi - pnd_L_hi0)**(p5) + pnd_L_c0 + &
             pnd_L_c1*(hi - pnd_L_hi0) + pnd_L_c2*(hi - pnd_L_hi0)**(c2)
      
      end function pond_log_L

!=======================================================================

      function pond_log_a0(hi) result(a0)

         ! Polynomial approximation of a0 for logistic pond hypsometry

         real(kind=dbl_kind), intent(in) :: &
            hi    ! ice thickness (m)
         
         real(kind=dbl_kind) :: &
            a0    ! location parameter for logistic hypsometry
         
         character(len=*),parameter :: subname='(pond_log_a0)'

         a0 = pnd_a0_c03*(hi - pnd_a0_hi0)**(p333) + &
              pnd_a0_c02*(hi - pnd_a0_hi0)**(p5) + pnd_a0_c0 + &
              pnd_a0_c1*(hi - pnd_a0_hi0)
      
      end function pond_log_a0

!=======================================================================

      function pond_log_y(hi) result(y)

         ! Polynomial approximation of y for logistic pond hypsometry

         real(kind=dbl_kind), intent(in) :: &
            hi    ! ice thickness (m)
         
         real(kind=dbl_kind) :: &
            y     ! vertical offset parameter for logistic hypsometry

         character(len=*),parameter :: subname='(pond_log_y)'

         y = pnd_y_c03*(hi - pnd_y_hi0)**(p333) + &
             pnd_y_c02*(hi - pnd_y_hi0)**(p5) + pnd_y_c0 + &
             pnd_y_c1*(hi - pnd_y_hi0)
      
      end function pond_log_y

!=======================================================================

      function pond_log_e(a, L, k, a0, y) result(e)

         ! Logistic hypsometric curve

         real(kind=dbl_kind), intent(in) :: &
            a, &  ! area fraction
            L, &  ! amplitude parameter for logistic hypsometry
            k, &  ! steepness parameter for logistic hypsometry
            a0, & ! location parameter for logistic hypsometry
            y     ! vertical offset parameter for logistic hypsometry
         
         real(kind=dbl_kind) :: &
            e     ! elevation above mean ice base on hypsometric curve at a (m)

         character(len=*),parameter :: subname='(pond_log_e)'

         e = L/(c1 + EXP(-k*(a-a0))) + y
      
      end function pond_log_e

!=======================================================================

      function pond_log_inte(a, L, k, a0, y) result(inte)

         ! Integral of logistic hypsometric curve

         real(kind=dbl_kind), intent(in) :: &
            a, &  ! area fraction
            L, &  ! amplitude parameter for logistic hypsometry
            k, &  ! steepness parameter for logistic hypsometry
            a0, & ! location parameter for logistic hypsometry
            y     ! vertical offset parameter for logistic hypsometry
         
         real(kind=dbl_kind) :: &
            inte  ! integral of hypsometric curve at a

         character(len=*),parameter :: subname='(pond_log_inte)'

         inte = L*LOG(EXP(-k*(a - a0)) + c1)/k + (L + y)*a &
                - L*LOG(EXP(k*a0) + c1)/k
      
      end function pond_log_inte
   
!=======================================================================

      function pond_log_Vp(a, L, k, a0, y) result(Vp)

         ! Integral of logistic hypsometric curve

         real(kind=dbl_kind), intent(in) :: &
            a, &  ! area fraction
            L, &  ! amplitude parameter for logistic hypsometry
            k, &  ! steepness parameter for logistic hypsometry
            a0, & ! location parameter for logistic hypsometry
            y     ! vertical offset parameter for logistic hypsometry
         
         real(kind=dbl_kind) :: &
            Vp    ! pond volume per unit category area (m)

         character(len=*),parameter :: subname='(pond_log_Vp)'

         Vp = a*pond_log_e(a,L,k,a0,y) - pond_log_inte(a,L,k,a0,y)
      
      end function pond_log_Vp

!=======================================================================

      function pond_log_deda(a, L, k, a0, y) result(deda)

         ! Derivative of logistic hypsometric curve

         real(kind=dbl_kind), intent(in) :: &
            a, &  ! area fraction
            L, &  ! amplitude parameter for logistic hypsometry
            k, &  ! steepness parameter for logistic hypsometry
            a0, & ! location parameter for logistic hypsometry
            y     ! vertical offset parameter for logistic hypsometry
         
         real(kind=dbl_kind) :: &
            deda  ! derivative of logistic hypsometric curve

         character(len=*),parameter :: subname='(pond_log_deda)'

         deda = k*(pond_log_e(a,L,k,a0,y) - y) &
                *(c1 - (pond_log_e(a,L,k,a0,y) - y)/L)
      
      end function pond_log_deda

!=======================================================================

! Compute the elevation of pond surface above the ice bottom

      subroutine pond_head(apondn, hpondn, alvln, hin, hpsurf)

      real (kind=dbl_kind), intent(in) :: &
         hin     , & ! category mean ice thickness
         apondn  , & ! pond area fraction of the category (incl. deformed)
         alvln   , & ! level ice fraction of category area
         hpondn      ! mean pond depth (m)

      real (kind=dbl_kind), intent(out) :: &
         hpsurf   ! height of pond surface above base of the ice (m)
      
      real (kind=dbl_kind) :: &
         L       , & ! amplitude parameter for logistic hypsometry
         a0      , & ! location parameter for logistic hypsometry
         y           ! vertical offset parameter for logistic hypsometry
      
      character(len=*),parameter :: subname='(pond_head)'

      if (trim(pndhead) == 'perched') then
         hpsurf = hin + hpondn
      elseif (trim(pndhead) == 'hyps') then
         if (trim(pndhyps) == 'fixed') then
            ! Applying a fixed aspect ratio to the ponds implicitly assumes
            ! that the hypsometric curve has a constant slope of double the
            ! aspect ratio. We'll assume that the deformed ice in the category
            ! also has the same mean thickness as the entire category.
            ! With these assumptions, we can derive the height of the mean
            ! pond surface above the mean base of the category
            if (apondn < alvln) then
               hpsurf = hin + c2*hpondn - alvln*pndaspect
            else
               hpsurf = hin + hpondn ! ponds cover all available area
            endif
         elseif (trim(pndhyps) == 'sealevel_lin') then
            ! Same as fixed aspect ratio but we have to calculate aspect ratio
            pndaspect = hin*(rhow - rhosi) / & 
                        (rhofresh*apond_sl**c2 - c2*rhow*apond_sl + rhow)
            if (apondn < alvln) then
               hpsurf = hin + c2*hpondn - alvln*pndaspect
            else
               hpsurf = hin + hpondn ! ponds cover all available area
            endif
         elseif (trim(pndhyps) == 'sealevel_log') then
            ! See above, logistic hypsometry
            if (hin < pnd_hi_min) then
               ! Same as fixed aspect ratio but we have to calculate aspect ratio
               pndaspect = hin*(rhow - rhosi) / & 
                           (rhofresh*apond_sl**c2 - c2*rhow*apond_sl + rhow)
               if (apondn < alvln) then
                  hpsurf = hin + c2*hpondn - alvln*pndaspect
               else
                  hpsurf = hin + hpondn ! ponds cover all available area
               endif
            else ! use logistic hypsomtery
               if (apondn < c1) then
                  L = pond_log_L(hin)
                  a0 = pond_log_a0(hin)
                  y = pond_log_y(hin)
                  hpsurf = pond_log_e(apondn,L,pnd_k,a0,y)
               else
                  hpsurf = hin + hpondn
               endif
            endif
         else
            call icepack_warnings_add(subname//" unsupported pndhyps option" )
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            if (icepack_warnings_aborted(subname)) return
         endif
      else
         call icepack_warnings_add(subname//" invalid pndhead option" )
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         if (icepack_warnings_aborted(subname)) return
      endif

      end subroutine pond_head

!=======================================================================

      end module icepack_meltpond_lvl

!=======================================================================
