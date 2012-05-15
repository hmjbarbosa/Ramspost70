!
!     *****************************************************************
      subroutine gradr(n1,n2,n3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,topt &
        ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim,ihtran)
      dimension vc3da(n1,n2,n3),vc3db(n1,n2,n3),e(1),topt(n1,n2), &
        xm(*),xt(*),ym(*),yt(*),zm(*),zt(*),dzm(*),dzt(*), &
        vctr1(*),vctr2(*)
      character*(*) dir,gpnt
      character*6 optyp
!
      optyp = 'gradnt'
      go to 10

      entry divcartr(n1,n2,n3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,topt &
        ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim)
      optyp = 'divcrt'
      go to 10   

      entry divstarr(n1,n2,n3,ia,iz,ja,jz,vc3da,vc3db,dir,gpnt,topt &
        ,xm,xt,ym,yt,zm,zt,deltay,dzm,dzt,vctr1,vctr2,ztop,jdim)
      optyp = 'divstr'

   10 continue

      jaa = ja
      jzz = jz
      if(jdim.eq.0) then
         jaa = 1
         jzz = 1
      endif

!      print*,'in grad-',dir,gpnt,optyp,'   !!!!!!!!!!!!!!!!'

      isizexy = n1 * n2
      isizeall = isizexy * 33
      iaddr=iralloc(isizeall,e,ioff)

      iglat   = ioff    + 1
      iglon   = iglat   + isizexy
      ifmapu  = iglon   + isizexy
      ifmapv  = ifmapu  + isizexy
      ifmapt  = ifmapv  + isizexy
      ifmapm  = ifmapt  + isizexy
      ifmapui = ifmapm  + isizexy
      ifmapvi = ifmapui + isizexy
      ifmapti = ifmapvi + isizexy
      ifmapmi = ifmapti + isizexy
      idxu    = ifmapmi + isizexy
      idxv    = idxu    + isizexy
      idxt    = idxv    + isizexy
      idxm    = idxt    + isizexy
      idyu    = idxm    + isizexy
      idyv    = idyu    + isizexy
      idyt    = idyv    + isizexy
      idym    = idyt    + isizexy
      itopu   = idym    + isizexy
      itopv   = itopu   + isizexy
      itopm   = itopv   + isizexy
      irtgt   = itopm   + isizexy
      irtgu   = irtgt   + isizexy
      irtgv   = irtgu   + isizexy
      irtgm   = irtgv   + isizexy
      if13u   = irtgm   + isizexy
      if13v   = if13u   + isizexy
      if13t   = if13v   + isizexy
      if13m   = if13t   + isizexy
      if23u   = if13m   + isizexy
      if23v   = if23u   + isizexy
      if23t   = if23v   + isizexy
      if23m   = if23t   + isizexy
      
      if(ihtran.eq.1)then
         CALL R_POLARST(N1,N2,E(IGLAT),E(IGLON),E(IFMAPU),E(IFMAPV) &
           ,E(IFMAPT),E(IFMAPM),E(IFMAPUI),E(IFMAPVI) &
           ,E(IFMAPTI),E(IFMAPMI),xm,xt,ym,yt,polelat,polelon,ihtran)
      else
         CALL AE0(N1*N2,E(IFMAPU),1.)
         CALL AE0(N1*N2,E(IFMAPV),1.)
         CALL AE0(N1*N2,E(IFMAPT),1.)
         CALL AE0(N1*N2,E(IFMAPM),1.)
         CALL AE0(N1*N2,E(IFMAPUI),1.)
         CALL AE0(N1*N2,E(IFMAPVI),1.)
         CALL AE0(N1*N2,E(IFMAPTI),1.)
         CALL AE0(N1*N2,E(IFMAPMI),1.)
      ENDIF

      CALL R_GRDSPC(N1,N2,E(IDXU),E(IDXV),E(IDXT),E(IDXM) &
        ,E(IDYU),E(IDYV),E(IDYT),E(IDYM) &
        ,E(IFMAPU),E(IFMAPV),E(IFMAPT),E(IFMAPM) &
        ,xm,xt,ym,yt,jdim,deltay)

      CALL R_TRANSFM(n1,n2,n3,TOPT,E(ITOPU),E(ITOPV),E(ITOPM) &
        ,E(IRTGT),E(IRTGU),E(IRTGV),E(IRTGM),E(IF13U),E(IF13V) &
        ,E(IF13T),E(IF13M),E(IF23U),E(IF23V),E(IF23T) &
        ,E(IF23M),E(IFMAPU),E(IFMAPV),E(IDXU),E(IDXV) &
        ,E(IDXT),E(IDXM),E(IDYU),E(IDYV),E(IDYT),E(IDYM) &
        ,xm,xt,ym,yt,zm,zt,jdim,ztop,ht,hw)

      if(dir.eq.'xdir')then
         if(gpnt.eq.'upnt')then
            call gradxur(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgu),e(irtgt),e(idxt),dzt,e(ifmapui),e(ifmapt) &
                ,e(if13t),hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'vpnt')then
            call gradxtr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgv),e(irtgm),e(idxm),dzt,e(ifmapvi),e(ifmapm) &
                ,e(if13m),hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'wpnt')then
            call gradxtr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgt),e(irtgu),e(idxu),dzm,e(ifmapti),e(ifmapu) &
                ,e(if13u),ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'tpnt')then
            call gradxtr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgt),e(irtgu),e(idxu),dzt,e(ifmapti),e(ifmapu) &
                ,e(if13u),hw,vctr2,'t',jdim)
         elseif(gpnt.eq.'npnt')then
            call gradxtr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgv),e(irtgm),e(idxm),dzm,e(ifmapvi),e(ifmapm) &
                ,e(if13m),ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'opnt')then
            call gradxur(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgu),e(irtgt),e(idxt),dzm,e(ifmapui),e(ifmapt) &
                ,e(if13t),ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'ppnt')then
            call gradxur(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgm),e(irtgv),e(idxv),dzt,e(ifmapmi),e(ifmapv) &
                ,e(if13v),hw,vctr2,'t',jdim)
        elseif(gpnt.eq.'mpnt')then
           call gradxur(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgm),e(irtgv),e(idxv),dzm,e(ifmapmi),e(ifmapv) &
                ,e(if13v),ht,vctr2,'w',jdim)
        endif
      elseif(dir.eq.'ydir')then
        if(gpnt.eq.'upnt')then
           call gradytr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgu),e(irtgm),e(idym),dzt,e(ifmapui),e(ifmapm) &
                ,e(if23m),hw,vctr2,'t',jdim)
        elseif(gpnt.eq.'vpnt')then
           call gradyvr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgv),e(irtgt),e(idyt),dzt,e(ifmapvi),e(ifmapt) &
                ,e(if23t),hw,vctr2,'t',jdim)
        elseif(gpnt.eq.'wpnt')then
           call gradytr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgt),e(irtgv),e(idyv),dzm,e(ifmapti),e(ifmapv) &
                ,e(if23v),ht,vctr2,'w',jdim)
        elseif(gpnt.eq.'tpnt')then
           call gradytr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgt),e(irtgv),e(idyv),dzt,e(ifmapti),e(ifmapv) &
                ,e(if23v),hw,vctr2,'t',jdim)
        elseif(gpnt.eq.'npnt')then
           call gradyvr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgv),e(irtgt),e(idyt),dzm,e(ifmapvi),e(ifmapt) &
                ,e(if23t),ht,vctr2,'w',jdim)
         elseif(gpnt.eq.'opnt')then
            call gradytr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgu),e(irtgm),e(idym),dzm,e(ifmapui),e(ifmapm) &
                ,e(if23m),ht,vctr2,'w',jdim)
        elseif(gpnt.eq.'ppnt')then
           call gradyvr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgm),e(irtgu),e(idyu),dzt,e(ifmapmi),e(ifmapu) &
                ,e(if23u),hw,vctr2,'t',jdim)
        elseif(gpnt.eq.'mpnt')then
           call gradyvr(n1,n2,n3,ia,iz,jaa,jzz,optyp,vc3da,vc3db,vctr1 &
                ,e(irtgm),e(irtgu),e(idyu),dzm,e(ifmapmi),e(ifmapu) &
                ,e(if23u),ht,vctr2,'w',jdim)
         endif
      elseif(dir.eq.'zdir')then
         if(gpnt.eq.'upnt')then
           call gradztr(n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgu),dzm)
         elseif(gpnt.eq.'vpnt')then
           call gradztr(n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgv),dzm)
         elseif(gpnt.eq.'wpnt')then
           call gradzwr(n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgt),dzt)
         elseif(gpnt.eq.'tpnt')then
           call gradztr(n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgt),dzm)
         elseif(gpnt.eq.'npnt')then
           call gradzwr(n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgv),dzt)
         elseif(gpnt.eq.'opnt')then
           call gradzwr(n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgu),dzt)
         elseif(gpnt.eq.'ppnt')then
           call gradztr(n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgm),dzm)
         elseif(gpnt.eq.'mpnt')then
           call gradzwr(n1,n2,n3,ia,iz,jaa,jzz,vc3da,vc3db,e(irtgm),dzt)
         endif
      endif

      ierr=irfree(iaddr)

      return
      end

!
!     *****************************************************************
!
!     ******************************************************************
!
      subroutine grad1(n1,n2,n3)
      dimension vc3da(n1,n2,n3),vc3db(n1,n2,n3),vc1da(*) &
       ,rtge(n1,n2),rtgc(n1,n2),dx(n1,n2),dy(n1,n2) &
       ,fmap(n1,n2),fmapi(n1,n2),dz(*),fq(n1,n2),hq(*),hq4(*)
      character*(*) optyp,lev
!
!     this is a general subroutine which computes any component of the
!     gradient or divergence of vc3da and stores it in vc3db.
!
      entry gradxur(n1,n2,n3,ia,iz,ja,jz,optyp,vc3da,vc3db,vc1da &
        ,rtge,rtgc,dx,dz,fmapi,fmap,fq,hq,hq4,lev,jd)

      if(optyp.eq.'gradnt')then
         do j = ja,jz
            do i = ia,iz
               do k = 1,n3
                  vc3db(i,j,k) = (vc3da(i,j,k) * rtge(i,j) &
                    - vc3da(i-1,j,k) * rtge(i-1,j)) &
                    * dx(i,j) / rtgc(i,j)
               enddo
            enddo
         enddo
      else
         do j = ja,jz
            do i = ia,iz
               do k = 1,n3
                  vc3db(i,j,k) = (vc3da(i,j,k) * rtge(i,j) * fmapi(i,j) & 
                    - vc3da(i-1,j,k) * rtge(i-1,j) * fmapi(i-1,j)) &
                    * dx(i,j) / rtgc(i,j) * fmap(i,j)
               enddo
            enddo
         enddo
      endif

      if(optyp.ne.'divstr')then
         if(lev.eq.'w')then
            do k = 1,n3
               hq4(k) = 0.25 * hq(k)
            enddo
         else
            do k = 2,n3
               hq4(k) = 0.25 * hq(k-1)
            enddo
         endif

         do j = ja,jz
            do i = ia,iz
               do k = 2,n3
                  vc1da(k) = hq4(k) * (vc3da(i,j,k) + vc3da(i,j,k-1) &
                    + vc3da(i-1,j,k) + vc3da(i-1,j,k-1))
              enddo
              do k = 2,n3-1
                 vc3db(i,j,k) = vc3db(i,j,k) &
                    + fq(i,j) * dz(k) * (vc1da(k+1) - vc1da(k))
               enddo
               vc3db(i,j,1) = vc3db(i,j,2)
               if(lev.eq.'w') vc3db(i,j,n3-1) = vc3db(i,j,n3-2)
               if(lev.eq.'t') vc3db(i,j,n3) = vc3db(i,j,n3-1)
            enddo
         enddo
      endif

      return
!
      entry gradxtr(n1,n2,n3,ia,iz,ja,jz,optyp,vc3da,vc3db,vc1da &
         ,rtge,rtgc,dx,dz,fmapi,fmap,fq,hq,hq4,lev,jd)

      if(optyp.eq.'gradnt')then
         do j = ja,jz
            do i = ia,iz
               do k = 1,n3
                  vc3db(i,j,k) = (vc3da(i+1,j,k) * rtge(i+1,j) &
                    - vc3da(i,j,k) * rtge(i,j)) &
                    * dx(i,j) / rtgc(i,j)
               enddo
            enddo
         enddo
      else
         do j = ja,jz
            do i = ia,iz
               do k = 1,n3
                  vc3db(i,j,k) = (vc3da(i+1,j,k) * rtge(i+1,j) &
                    * fmapi(i+1,j) - vc3da(i,j,k) * rtge(i,j) &
                    * fmapi(i,j)) * dx(i,j) / rtgc(i,j) * fmap(i,j)
               enddo
            enddo
         enddo
      endif

      if(optyp.ne.'divstr')then
         if(lev.eq.'w')then
            do k = 1,n3
               hq4(k) = 0.25 * hq(k)
            enddo
         else
            do k = 2,n3
               hq4(k) = 0.25 * hq(k-1)
            enddo
         endif

         do j = ja,jz
            do i = ia,iz
               do k = 2,n3
                  vc1da(k) = hq4(k) * (vc3da(i,j,k) + vc3da(i,j,k-1) &
                    + vc3da(i+1,j,k) + vc3da(i+1,j,k-1))
              enddo
              do k = 2,n3-1
                 vc3db(i,j,k) = vc3db(i,j,k) &
                    + fq(i,j) * dz(k) * (vc1da(k+1) - vc1da(k))
               enddo
               vc3db(i,j,1) = vc3db(i,j,2)
               if(lev.eq.'w') vc3db(i,j,n3-1) = vc3db(i,j,n3-2)
               if(lev.eq.'t') vc3db(i,j,n3) = vc3db(i,j,n3-1)
            enddo
         enddo
      endif

      return
!
      entry gradyvr(n1,n2,n3,ia,iz,ja,jz,optyp,vc3da,vc3db,vc1da &
          ,rtge,rtgc,dy,dz,fmapi,fmap,fq,hq,hq4,lev,jd)

      if(optyp.eq.'gradnt')then
        do j = ja,jz
           do i = ia,iz
              do k = 1,n3
                 vc3db(i,j,k) = (vc3da(i,j,k) * rtge(i,j) &
                    - vc3da(i,j-jd,k) * rtge(i,j-jd)) &
                    * dy(i,j) / rtgc(i,j)
               enddo
            enddo
         enddo
      else
         do j = ja,jz
            do i = ia,iz
               do k = 1,n3
                  vc3db(i,j,k) = (vc3da(i,j,k) * rtge(i,j) * fmapi(i,j) &
                    - vc3da(i,j-jd,k) * rtge(i,j-jd) * fmapi(i,j-jd)) &
                    * dy(i,j) / rtgc(i,j) * fmap(i,j)
               enddo
            enddo
         enddo
      endif

      if(optyp.ne.'divstr')then
         if(lev.eq.'w')then
            do k = 1,n3
               hq4(k) = 0.25 * hq(k)
            enddo
         else
            do k = 2,n3
               hq4(k) = 0.25 * hq(k-1)
            enddo
         endif

         do j = ja,jz
            do i = ia,iz
               do k = 2,n3
                  vc1da(k) = hq4(k) * (vc3da(i,j,k) + vc3da(i,j,k-1) &
                    + vc3da(i,j-jd,k) + vc3da(i,j-jd,k-1))
              enddo
              do k = 2,n3-1
                 vc3db(i,j,k) = vc3db(i,j,k) &
                    + fq(i,j) * dz(k) * (vc1da(k+1) - vc1da(k))
               enddo
               vc3db(i,j,1) = vc3db(i,j,2)
               if(lev.eq.'w') vc3db(i,j,n3-1) = vc3db(i,j,n3-2)
               if(lev.eq.'t') vc3db(i,j,n3) = vc3db(i,j,n3-1)
            enddo
         enddo
      endif

      return
!
      entry gradytr(n1,n2,n3,ia,iz,ja,jz,optyp,vc3da,vc3db,vc1da &
                 ,rtge,rtgc,dy,dz,fmapi,fmap,fq,hq,hq4,lev,jd)

      if(optyp.eq.'gradnt')then
         do j = ja,jz
            do i = ia,iz
               do k = 1,n3
                  vc3db(i,j,k) = (vc3da(i,j+jd,k) * rtge(i,j+jd) &
                    - vc3da(i,j,k) * rtge(i,j)) &
                    * dy(i,j) / rtgc(i,j)
              enddo
           enddo
        enddo
      else
        do j = ja,jz
           do i = ia,iz
              do k = 1,n3
                 vc3db(i,j,k) = (vc3da(i,j+jd,k) * rtge(i,j+jd) &
                    * fmapi(i,j+jd) &
                    - vc3da(i,j,k) * rtge(i,j) * fmapi(i,j)) &
                    * dy(i,j) / rtgc(i,j) * fmap(i,j)
               enddo
            enddo
         enddo
      endif

      if(optyp.ne.'divstr')then
         if(lev.eq.'w')then
            do k = 1,n3
               hq4(k) = 0.25 * hq(k)
            enddo
         else
            do k = 2,n3
               hq4(k) = 0.25 * hq(k-1)
            enddo
         endif

         do j = ja,jz
            do i = ia,iz
               do k = 2,n3
                  vc1da(k) = hq4(k) * (vc3da(i,j,k) + vc3da(i,j,k-1) &
                    + vc3da(i,j+jd,k) + vc3da(i,j+jd,k-1))
              enddo
              do k = 2,n3-1
                 vc3db(i,j,k) = vc3db(i,j,k) &
                    + fq(i,j) * dz(k) * (vc1da(k+1) - vc1da(k))
               enddo
               vc3db(i,j,1) = vc3db(i,j,2)
               if(lev.eq.'w') vc3db(i,j,n3-1) = vc3db(i,j,n3-2)
               if(lev.eq.'t') vc3db(i,j,n3) = vc3db(i,j,n3-1)
            enddo
         enddo
      endif

      return
!
      entry gradzwr(n1,n2,n3,ia,iz,ja,jz,vc3da,vc3db,rtgc,dz)

      do j = ja,jz
         do i = ia,iz
            do k = 2,n3
               vc3db(i,j,k) = (vc3da(i,j,k) - vc3da(i,j,k-1)) * dz(k) &
                 / rtgc(i,j)
            enddo
         enddo
      enddo
      return
!     
      entry gradztr(n1,n2,n3,ia,iz,ja,jz,vc3da,vc3db,rtgc,dz)

      do j = ja,jz
         do i = ia,iz
            do k = 1,n3-1
               vc3db(i,j,k) = (vc3da(i,j,k+1) - vc3da(i,j,k)) * dz(k) &
                 / rtgc(i,j)
            enddo
         enddo
      enddo
      return
      end
!     ******************************************************************
!
      SUBROUTINE XYTOPS(X,Y,PLA,PLO,ERAD)
!
!     This convert x,y-polar stereographic coordinates to
!     lat/lon values.
!     longitude:   0 - 360  ; positive to the east
!     latitude : -90 -  90  ; positive for northern hemisphere
!     It is assumed that the x-axis point towards the east why the
!     longitude is rotated relative to the 'standard pol.ste.' location
!            

      PI180=3.14159/180.
      VDIST = ERAD*2.0
!
!     calculate distance from (0,0) and define (0,0) as 90,0 (90,-90 in
!     the rotated system)
!
      DIST=(X**2+Y**2)**0.5

      IF(DIST.EQ.0) THEN
         PLA= 90.0
         PLO=-90.0
      ELSE
!
!     calculate the latitude by means of atan
!
         PLA=ATAN(DIST/VDIST)/PI180
         PLA=90.0-2.0*PLA
!
!     calculate the longitude taking the directions into account
!
         IF(X.EQ.0.0) THEN
            IF(Y.GT.0.0) THEN
               PLO= 90.0
            ELSE
               PLO=-90.0
            END IF
         ELSE
            IF(X.GT.0.0) THEN
               PLO=ATAN(Y/X)/PI180
            ELSE
               PLO=ATAN(Y/X)/PI180+180.0
            END IF
         END IF
      END IF
!
!     rotate the longitude
!
      PLO=AMOD(PLO+450.0,360.0)
      RETURN
      END
!-------------------------------------------------------------------
!      subroutine ll_xy(xlat,xlon,plat,plon,x,y)
!      call getops(pla,plo,xlat,xlon,plat,plon)
!      call pstoxy(x,y,pla,plo,6376000.)
!      return
!      end
!
!
!      subroutine xy_ll(xlat,xlon,plat,plon,x,y)
!      call xytops(x,y,pla,plo,6376000.)
!      call pstoge(pla,plo,xlat,xlon,plat,plon)
!      return
!      end
!
!      function hg_xy(nih,xh,hg)
!      dimension xh(*)
!      ix=int(hg)
!      hg_xy=xh(ix)+(xh(ix+1)-xh(ix))*(hg-ix)
!      return
!      end
!
!      function hifromx(nih,xh,x)
!      dimension xh(*)
!
!!      print*,'hifromx-',x,nih,xh(1),xh(nih)
!      if(x.lt.xh(1).or.x.gt.xh(nih))then
!         print*, 'x, y, or z value exceeds hgrid limits'
!         stop 'xtohi'
!      endif
!
!      ilow=1
!      ihigh=nih
!      do m=2,11
!         imid=(ilow+ihigh)/2
!         msw=int(max(0.,min(1.5,1.e20*(x-xh(imid)))))
!         ilow=ilow*(1-msw)+imid*msw
!         ihigh=ihigh*msw+imid*(1-msw)
!      enddo
!      hifromx=float(ilow)+(x-xh(ilow))/(xh(ihigh)-xh(ilow))
!      end


!
!      function hg_z(nih,njh,nkh,xh,yh,zh,hgx,hgy,hgz,topth,nh1,nh2,ztop
!     +     ,topo)
!      dimension topth(nh1,nh2),xh(*),yh(*),zh(*)
!
!      ihp=int(hgx)
!      jhp=int(hgy)
!      qh20=hgx-float(ihp)
!      qh02=hgy-float(jhp)
!      topo=(1.0-qh02)*((1.0-qh20)*topth(ihp  ,jhp  )
!     +     +qh20 *topth(ihp+1,jhp  ))
!     +     +qh02 *((1.0-qh20)*topth(ihp  ,jhp+1)
!     +     +qh20 *topth(ihp+1,jhp+1))
!      rtg=1.-topo/ztop
!
!
!      iz=int(hgz)
!      hg_z=topo+(zh(iz)+(zh(iz+1)-zh(iz))*(hgz-iz))*rtg
!      print*,'hg_z-',hgx,hgy,hgz,topo,hg_z,topth(ihp,jhp)
!      return
!      end
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE UVTOUV(UP,VP,UE,VE,X,Y,DIR,RLAT,WLON1,ERAD)
!
!     Converts u,v in polar stereographic space to u,v in normal space
!     for DIR='STOG' else the convertion will be from normal to polar
!     stereographic space.
!     The conversion is made by transforming to dd,ff and rotate dd
!     after which the new dd,ff is converted back to u,v
!     NOTE: dd is not the normal dd but the angle between the x-axis
!           and the wind vector
!
!     TSP 03 JULY 89
!
      CHARACTER*4 DIR
      PI180=3.14159/180.0
!
!     Determine the location in pla,plo- and gla,glo coordinates
!
      CALL XYTOPS(X,Y,PLA,PLO,ERAD)
      CALL PSTOGE(PLA,PLO,GLA,GLO,RLAT,WLON1)
!
!     Determine the rotaion angle alpha which depends on the sign of RLAT
!
      HSIGN=1.0
      IF(RLAT.LT.0.0) HSIGN=-1.0
      ARGA1=COS(PLA*PI180)*COS(HSIGN*GLA*PI180)
      IF(ARGA1.NE.0.0.AND.((X**2+Y**2).GT.1.0)) THEN
         ARGA2=(SIN(HSIGN*RLAT*PI180)- &
               SIN(PLA*PI180)*SIN(HSIGN*GLA*PI180))/ARGA1
      ELSE
         ARGA2=1.0
      END IF
      IF(ARGA2.GT.1.0) THEN
         ARGA2=1.0
      ELSE IF(ARGA2.LT.-1.0) THEN
         ARGA2=-1.0
      END IF
      ALPHA=ACOS(ARGA2)
!
!     Determine rotation angle and sin/cos values between pol.ste. and
!     the geo. system
!
      ROTANG=(PLO+270.0)*PI180
      RSIN = SIN(ROTANG)
      RCOS = COS(ROTANG)
!
      IF(DIR.EQ.'STOG') THEN
!
!     Convert up,vp to us,vs (ie from x,y directions in pol.ste. to
!     east, north in the pla, plo pol.ste. rotated system
!
        US   =-RSIN*UP + RCOS*VP
        VS   =-RCOS*UP - RSIN*VP
!
!     Convert us,vs to DD,FF with DD <-180, 180>
!         
         FF=SQRT(US**2+VS**2)
         IF(US.EQ.0.0.AND.VS.GE.0.0) THEN
            DD=PI180*90.0
         ELSE IF(US.EQ.0.0.AND.VS.LT.0.0) THEN
            DD=PI180*270.0
         ELSE IF(US.GT.0.0) THEN
            DD=ATAN(VS/US)
         ELSE
            DD=ATAN(VS/US)
            DD=DD+PI180*180.0
         END IF
!
!     The actual rotation depends on RLAT and PLO
!
         IF(RLAT.GE.0.0) THEN
            IF(PLO.LT.180.0) THEN
               DD=DD+ALPHA
            ELSE
               DD=DD-ALPHA
            END IF
         ELSE
            IF(PLO.LT.180.0) THEN
               DD=DD-ALPHA+PI180*180.0
            ELSE
               DD=DD+ALPHA+PI180*180.0
            END IF
         END IF
!
!     Conversion from dd,ff to ue,ve
!
         UE=FF*COS(DD)
         VE=FF*SIN(DD)
!
!     This corresponds to the conversion from normal to polar stereographic
!         
      ELSE
!         
!     Convert UE,VE to DD,FF with DD <-180, 180>
!         
         FF=SQRT(UE**2+VE**2)
         IF(UE.EQ.0.0.AND.VE.GE.0.0) THEN
            DD=PI180*90.0
         ELSE IF(UE.EQ.0.0.AND.VE.LT.0.0) THEN
            DD=PI180*270.0
         ELSE IF(UE.GT.0.0) THEN
            DD=ATAN(VE/UE)
         ELSE
            DD=ATAN(VE/UE)
            DD=DD+PI180*180.0
         END IF
!
!     The actual rotation depends on RLAT and PLO
!
         IF(RLAT.GE.0.0) THEN
            IF(PLO.LT.180.0) THEN
               DD=DD-ALPHA
            ELSE
               DD=DD+ALPHA
            END IF
         ELSE
            IF(PLO.LT.180.0) THEN
               DD=DD+ALPHA-PI180*180.0
            ELSE
               DD=DD-ALPHA-PI180*180.0
            END IF
         END IF
!
!     Conversion from dd,ff to us,vs
!
         US=FF*COS(DD)
         VS=FF*SIN(DD)
!
!     Convert us,vs to up,vp (ie from east, north directions in pol.ste. to
!     x,y in the pol.ste. rotated system
!
        UP   =-RSIN*US - RCOS*VS
        VP   = RCOS*US - RSIN*VS
      END IF
!
      RETURN
      END
!
!     ******************************************************************
!
      SUBROUTINE R_TRANSFM(N1,N2,n3,TOPT,TOPU,TOPV,TOPM,RTGT,RTGU,RTGV &
        ,RTGM,F13U,F13V,F13T,F13M,F23U,F23V,F23T,F23M &
        ,FMAPU,FMAPV,DXU,DXV,DXT,DXM,DYU,DYV,DYT,DYM &
        ,xm,xt,ym,yt,zm,zt,jdim,ztop,ht,hw)
      DIMENSION TOPT(N1,N2),TOPU(N1,N2),TOPV(N1,N2),TOPM(N1,N2) &
        ,RTGT(N1,N2),RTGU(N1,N2),RTGV(N1,N2),RTGM(N1,N2) &
        ,F13U(N1,N2),F13V(N1,N2),F13T(N1,N2),F13M(N1,N2) &
        ,F23U(N1,N2),F23V(N1,N2),F23T(N1,N2),F23M(N1,N2) &
        ,FMAPU(N1,N2),FMAPV(N1,N2) &
        ,DXU(N1,N2),DXV(N1,N2),DXT(N1,N2),DXM(N1,N2) &
        ,DYU(N1,N2),DYV(N1,N2),DYT(N1,N2),DYM(N1,N2) &
        ,xm(*),xt(*),ym(*),yt(*),zm(*),zt(*),ht(*),hw(*)
      DATA TERDEV/0./
      SAVE
!
!    This routine computes the coordinate transformation constants
!    based on the topographical values of TOPT.
!
      DO J=1,N2
        DO I=1,N1-1
           TOPU(I,J)=TOPT(I,J)+(TOPT(I+1,J)-TOPT(I,J)) &
                    *(XM(I)-XT(I))/(XT(I+1)-XT(I))
           TERDEV=MAX(TERDEV,ABS(TOPT(I,J)))
        ENDDO
        TOPU(N1,J)=TOPT(N1,J)+(TOPT(N1,J)-TOPT(N1-1,J)) &
                   *(XM(N1)-XT(N1))/(XT(N1)-XT(N1-1))
      ENDDO
      IF(TERDEV.LT.1.E-6)THEN
         ITOPO=0
      ELSE
         ITOPO=1
      ENDIF
!
      IF(JDIM.EQ.1)THEN
         DO I=1,N1
            DO J=1,N2-1
               TOPV(I,J)=TOPT(I,J)+(TOPT(I,J+1)-TOPT(I,J)) &
                       *(YM(J)-YT(J))/(YT(J+1)-YT(J))
              TOPM(I,J)=TOPU(I,J)+(TOPU(I,J+1)-TOPU(I,J)) &
                       *(YM(J)-YT(J))/(YT(J+1)-YT(J))
           ENDDO
           TOPV(I,N2)=TOPT(I,N2)+(TOPT(I,N2)-TOPT(I,N2-1)) &
                   *(YM(N2)-YT(N2))/(YT(N2)-YT(N2-1))
           TOPM(I,N2)=TOPU(I,N2)+(TOPU(I,N2)-TOPU(I,N2-1)) &
                   *(YM(N2)-YT(N2))/(YT(N2)-YT(N2-1))
         ENDDO
      ELSE
         DO I=1,N1
            DO J=1,N2
               TOPV(I,J)=TOPT(I,J)
               TOPM(I,J)=TOPU(I,J)
            ENDDO
         ENDDO
      ENDIF
!
      DO K=1,N3
        HT(K)=ZT(K)/ZTOP-1.0
        HW(K)=ZM(K)/ZTOP-1.0
      ENDDO
!
      DO J=1,N2
        DO I=1,N1
          RTGT(I,J)=1.-TOPT(I,J)/ZTOP
          RTGU(I,J)=1.-TOPU(I,J)/ZTOP
          RTGV(I,J)=1.-TOPV(I,J)/ZTOP
          RTGM(I,J)=1.-TOPM(I,J)/ZTOP
          F13T(I,J)=(TOPU(I,J)-TOPU(I-1,J))*DXT(I,J)/RTGT(I,J)
          F13U(I,J)=(TOPT(I+1,J)-TOPT(I,J))*DXU(I,J)/RTGU(I,J)
          F13V(I,J)=(TOPM(I,J)-TOPM(I-1,J))*DXV(I,J)/RTGV(I,J)
          F13M(I,J)=(TOPV(I+1,J)-TOPV(I,J))*DXM(I,J)/RTGM(I,J)
          F23T(I,J)=(TOPV(I,J)-TOPV(I,J-JDIM))*DYT(I,J)/RTGT(I,J)
          F23U(I,J)=(TOPM(I,J)-TOPM(I,J-JDIM))*DYU(I,J)/RTGU(I,J)
          F23V(I,J)=(TOPT(I,J+JDIM)-TOPT(I,J))*DYV(I,J)/RTGV(I,J)
          F23M(I,J)=(TOPU(I,J+JDIM)-TOPU(I,J))*DYM(I,J)/RTGM(I,J)
        ENDDO
      ENDDO
      DO J=1,N2
        F13T(1,J)=F13U(1,J)
        F13V(1,J)=F13M(1,J)
        F13U(N1,J)=F13T(N1,J)
        F13M(N1,J)=F13V(N1,J)
      ENDDO
      IF(JDIM.EQ.1)THEN
        DO I=1,N1
          F23T(I,1)=F23V(I,1)
          F23U(I,1)=F23M(I,1)
          F23V(I,N2)=F23T(I,N2)
          F23M(I,N2)=F23U(I,N2)
        ENDDO
      ENDIF
!
      RETURN
      END
!     *****************************************************************
!
      SUBROUTINE R_POLARST(N1,N2,GLAT,GLON,FMAPU,FMAPV,FMAPT,FMAPM &
       ,FMAPUI,FMAPVI,FMAPTI,FMAPMI,xm,xt,ym,yt,polelat,polelon,ihtran)
      parameter(erad=6367000.)
      DIMENSION GLAT(N1,N2),GLON(N1,N2),FMAPU(N1,N2),FMAPUI(N1,N2) &
              ,FMAPV(N1,N2),FMAPT(N1,N2),FMAPM(N1,N2),FMAPVI(N1,N2) &
              ,FMAPTI(N1,N2),FMAPMI(N1,N2),xm(*),xt(*),ym(*),yt(*)

!     Calculates map factors at u,v,t,m-points and geographical lat/lon
!     for the t-points for a given polar stereographic grid

      c1 = (2. * erad) ** 2
      DO I=1,n1
         DO J=1,n2
            xm2 = xm(i) * xm(i)
            xt2 = xt(i) * xt(i)
            ym2 = ym(j) * ym(j)
            yt2 = yt(j) * yt(j)

            FMAPT(I,J) = 1. + (xt2 + yt2) / c1
            FMAPU(I,J) = 1. + (xm2 + yt2) / c1
            FMAPV(I,J) = 1. + (xt2 + ym2) / c1
            FMAPM(I,J) = 1. + (xm2 + ym2) / c1

            FMAPUI(I,J) = 1.0 / FMAPU(I,J)
            FMAPVI(I,J) = 1.0 / FMAPV(I,J)
            FMAPTI(I,J) = 1.0 / FMAPT(I,J)
            FMAPMI(I,J) = 1.0 / FMAPM(I,J)

            call xy_ll(GLAT(I,J),GLON(I,J),polelat,polelon,XT(I),YT(J))
         ENDDO
      ENDDO

      IF (IHTRAN .EQ. 0) THEN
         CALL AE0(N1*N2,FMAPU,1.)
         CALL AE0(N1*N2,FMAPV,1.)
         CALL AE0(N1*N2,FMAPT,1.)
         CALL AE0(N1*N2,FMAPM,1.)
         CALL AE0(N1*N2,FMAPUI,1.)
         CALL AE0(N1*N2,FMAPVI,1.)
         CALL AE0(N1*N2,FMAPTI,1.)
         CALL AE0(N1*N2,FMAPMI,1.)
      ENDIF
      RETURN
      END
!     *****************************************************************
!
      SUBROUTINE R_GRDSPC(N1,N2,DXU,DXV,DXT,DXM,DYU,DYV,DYT,DYM &
        ,FMAPU,FMAPV,FMAPT,FMAPM,xm,xt,ym,yt,jdim,deltay)
      DIMENSION DXU(N1,N2),DXV(N1,N2),DXT(N1,N2),DXM(N1,N2) &
        ,DYU(N1,N2),DYV(N1,N2),DYT(N1,N2),DYM(N1,N2) &
        ,FMAPU(N1,N2),FMAPV(N1,N2),FMAPT(N1,N2),FMAPM(N1,N2) &
        ,xm(*),xt(*),ym(*),yt(*)
!
      DO J=1,N2
         DO I=1,N1-1
            DXU(I,J)=FMAPU(I,J)/(XT(I+1)-XT(I))
            DXM(I,J)=FMAPM(I,J)/(XT(I+1)-XT(I))
         ENDDO
         DXU(N2,J)=DXU(N2-1,J)*FMAPU(N2,J)/FMAPU(N2-1,J)
         DXM(N2,J)=DXM(N2-1,J)*FMAPM(N2,J)/FMAPM(N2-1,J)       
         DO I=2,N1
            DXV(I,J)=FMAPV(I,J)/(XM(I)-XM(I-1))
            DXT(I,J)=FMAPT(I,J)/(XM(I)-XM(I-1))
         ENDDO
         DXV(1,J)=DXV(2,J)*FMAPV(1,J)/FMAPV(2,J)
         DXT(1,J)=DXT(2,J)*FMAPT(1,J)/FMAPT(2,J)       
      ENDDO
!
      IF(JDIM.EQ.1)THEN
         DO I=1,N1
            DO J=1,N2-1
               DYV(I,J)=FMAPV(I,J)/(YT(J+1)-YT(J))
               DYM(I,J)=FMAPM(I,J)/(YT(J+1)-YT(J))
            ENDDO
            DYV(I,N2)=DYV(I,N2-1)*FMAPV(I,N2)/FMAPV(I,N2-1)
            DYM(I,N2)=DYM(I,N2-1)*FMAPM(I,N2)/FMAPM(I,N2-1)       
            DO J=2,N2
               DYU(I,J)=FMAPU(I,J)/(YM(J)-YM(J-1))
               DYT(I,J)=FMAPT(I,J)/(YM(J)-YM(J-1))
            ENDDO
            DYU(I,1)=DYU(I,2)*FMAPU(I,1)/FMAPU(I,2)
            DYT(I,1)=DYT(I,2)*FMAPT(I,1)/FMAPT(I,2)       
         ENDDO
      ELSE
         DO I=1,N1
            DO J=1,N2
               DYU(I,J)=1./DELTAY
               DYV(I,J)=1./DELTAY
               DYT(I,J)=1./DELTAY
               DYM(I,J)=1./DELTAY
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END
!     ******************************************************************
!
      SUBROUTINE PSTOXY(X,Y,PLA,PLO,ERAD)
!
!     This program convert polar stereographic coordinates to x,y ditto
!     longitude:   0 - 360  ; positive to the east
!     latitude : -90 -  90  ; positive for northern hemisphere
!     it is assumed that the x-axis point towards the east and
!     corresponds to longitude = 0
!
!     TSP 20/06-89
!
!     constants and functions
!            
      FAC(PLA) = ERAD*2.0/(1.0+SIN(PLA*PI180))*COS(PLA*PI180)
      XC(PLA,PLO) = FAC(PLA)*COS(PLO*PI180)
      YC(PLA,PLO) = FAC(PLA)*SIN(PLO*PI180)      
      PI180=3.14159/180.0
!
!     Calculate the coordinates
!
      X = XC(PLA,PLO)
      Y = YC(PLA,PLO)
!
      RETURN
      END
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PSTOGE(PLA,PLO,GLAT,GLON,RLAT,WLON1)
!
!     Convert polar stereographic coordinates to geographical lat/lon
!     ditto with the pol.ste. pole at rlat,wlon1 (these names are
!     used to be compatible to the ramsin-parameters)
!     longitude:   0 ; 360 positive east (on input)
!               -180 ; 180 positive east (on output)
!     latitude : -90 ;  90 posive on northern hemisphere
!     It is assumed that the polar stereographic coordinates have been
!     rotated to the standard format with 0 degrees longitude along wlon1
!
!     TSP 21 JUNE 89
!
!     set flag for n/s hemisphere
!
      PI180=3.14159/180.
!      
      IF(RLAT.GE.0.0) THEN
         HSIGN= 1.0
      ELSE
         HSIGN=-1.0
      END IF
!
!     test for a n/s pole case
!
      IF(RLAT.EQ.90.0) THEN
	 GLAT=PLA
         GLON=MOD(PLO+WLON1,360.0)
         GO TO 2000
      END IF
      IF(RLAT.EQ.-90.0) THEN
         GLAT=-PLA
         GLON=MOD(PLO+WLON1,360.0)
         GO TO 2000
      END IF
!
!     test for longitude on 'greenwich or date line'
!
      IF(PLO.EQ.0) THEN
         GLAT=RLAT-90.0+PLA
         IF(GLAT.LT.-90.0) THEN
            GLAT=-180.0-GLAT
            GLON=MOD(WLON1+180.0,360.0)
         ELSE
            GLON=WLON1
         END IF
         GO TO 2000
      END IF      
      IF(PLO.EQ.180.0) THEN
         GLAT=RLAT+90.0-PLA
         IF(GLAT.GT.90.0) THEN
            GLAT=180.0-GLAT
            GLON=MOD(WLON1+180.0,360.0)
         ELSE
            GLON=WLON1
         END IF
         GO TO 2000         
      END IF
!
!     Determine longitude distance relative to wlon1 so it belongs to
!     the absolute interval 0 - 180
!
      ARGU1=PLO
      IF(PLO.GT.180.0) ARGU1 = PLO-360.0
!
!     Get the latitude, the help circle BB and the longitude by first
!     calculating the argument and legalize it - then take the inverse fct.
!
      IF(HSIGN.GT.0.0) THEN
         ARG2A = SIN(PLA*PI180)*SIN(HSIGN*RLAT*PI180)+ &
             COS(PLA*PI180)*COS(RLAT*PI180)*COS((180.0-ARGU1)*PI180)
      ELSE
        ARG2A = SIN(PLA*PI180)*SIN(HSIGN*RLAT*PI180)+ &
             COS(PLA*PI180)*COS(RLAT*PI180)*COS(ARGU1*PI180)
      END IF
      ARG2A = MIN(ARG2A, 1.0)
      ARG2A = MAX(ARG2A,-1.0)
      GLAT  = HSIGN*ASIN(ARG2A)
!
      IF(HSIGN.GT.0.0) THEN
         ARG2A = COS(RLAT*PI180)*SIN(PLA*PI180)+ &
             SIN(RLAT*PI180)*COS(PLA*PI180)*COS(ARGU1*PI180)
      ELSE
        ARG2A = COS(RLAT*PI180)*SIN(PLA*PI180)+ &
            SIN(-RLAT*PI180)*COS(PLA*PI180)*COS((180.0-ARGU1)*PI180)
      END IF
      ARG2A = MIN(ARG2A, 1.0)
      ARG2A = MAX(ARG2A,-1.0)      
      BB    = ACOS(ARG2A)
!
      ARG2A = COS(GLAT)*COS(BB)/(1.0-SIN(GLAT)**2)
      ARG2A = MIN(ARG2A, 1.0)
      ARG2A = MAX(ARG2A,-1.0)      
      GLON  = ACOS(ARG2A)
!     
!     convert the radians to degrees 
!
        GLAT = GLAT/PI180
        GLON = GLON/PI180
!
!       the solution is symmetric so the direction must be if'ed
!
        IF(ARGU1.LT.0.0) THEN
           GLON = 360.0-GLON
        END IF
        GLON=AMOD(GLON+WLON1,360.0)
!
 2000 CONTINUE
!
!     the resultant longitude must be in the interval from -180, 180
!      
      IF(GLON.GT.180.0) GLON=GLON-360.0
      RETURN
      END
!
!     ******************************************************************
!
!     ****************************************************************
!
      SUBROUTINE GETOPS(PLA,PLO,GLAT,GLON,RLAT,WLON1)
!
!     Convert geographical lat/lon coordinates to polar stereographic
!     ditto with the pol.ste. pole at RLAT,WLON1 (these names are
!     used to be compatible to the RAMSIN-parameters)
!     longitude:-180 ; 180 positive east (on input)
!              :   0 ; 360 positive east (on output)
!     latitude : -90 ;  90 posive on northern hemisphere
!     The result is rotated 270 degrees relative to 'standard pol.ste.'
!     WLON1 is defined in the same way as the input
!     approach so as to get the x-axis to point towards the east, and the
!     y-axis towards the north along 0 degrees (at NP south along 180)
!
!     TSP 20/06-89
      double precision pi180,c1,c2,c3,c4,c5,c6,c7,arg2a,bb,pla1,alpha &
        ,plo1,pla90,argu2
!
!     constants
!
      c1=1.
      PI180 = dasin(c1)/90.
!
!     Set flag for N/S hemisphere and convert longitude to <0 ; 360> interval
!
      IF(RLAT.GE.0.0) THEN
         HSIGN= 1.0
      ELSE
         HSIGN=-1.0
      END IF
      GLOR=GLON
      IF(GLOR.LT.0.0) GLOR=360.0+GLOR
      RWLON1=WLON1
      IF(RWLON1.LT.0.0) RWLON1=360.0+WLON1
!
!     Test for a N/S pole case
!
      IF(RLAT.EQ.90.0) THEN
         PLA=GLAT
         PLO=AMOD(GLOR+270.0,360.0)
         GO TO 2000
      END IF
      IF(RLAT.EQ.-90.0) THEN
         PLA=-GLAT
         PLO=AMOD(GLOR+270.0,360.0)
         GO TO 2000
      END IF
!
!     Test for longitude on 'Greenwich or date line'
!
      IF(GLOR.EQ.RWLON1) THEN
         IF(GLAT.GT.RLAT) THEN
            PLA=90.0-GLAT+RLAT
            PLO=90.0
         ELSE
            PLA=90.0-RLAT+GLAT
            PLO=270.0
         END IF
         GO TO 2000
      END IF      
      IF(AMOD(GLOR+180.0,360.0).EQ.RWLON1) THEN
         PLA=RLAT-90.0+GLAT
         IF(PLA.LT.-90.0) THEN
            PLA=-180.0-PLA
            PLO=270.0
         ELSE
            PLO= 90.0
         END IF
         GO TO 2000         
      END IF
!
!     Determine longitude distance relative to RWLON1 so it belongs to
!     the absolute interval 0 - 180
!
      ARGU1 = GLOR-RWLON1
      IF(ARGU1.GT. 180.0) ARGU1 = ARGU1-360.0
      IF(ARGU1.LT.-180.0) ARGU1 = ARGU1+360.0
!
!     1. Get the help circle BB and angle ALPHA (legalize arguments)
!
      c2=glat*pi180
      c3=argu1*pi180
      ARG2A = dCOS(c2)*dCOS(c3)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)         
      BB    = dACOS(ARG2A)
!
      c4=hsign*glat*pi180
      ARG2A = dSIN(c4)/dSIN(BB)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)
      ALPHA = dASIN(ARG2A)
!
!     2. Get PLA and PLO (still legalizing arguments)
!
      c5=rlat*pi180
      c6=hsign*rlat*pi180
      ARG2A = dCOS(c5)*dCOS(BB)+ &
             dSIN(c6)*dSIN(c4)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)         
      PLA1   = dASIN(ARG2A)
!
      ARG2A = dSIN(BB)*dCOS(ALPHA)/dCOS(PLA1)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)
      PLO1   = dASIN(ARG2A)
!
!    Test for passage of the 90 degree longitude (duallity in PLO)
!         Get PLA for which PLO=90 when GLAT is the latitude
!
      ARG2A = dSIN(c4)/dSIN(c6)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)         
      PLA90 = dASIN(ARG2A)
!
!         Get help arc BB and angle ALPHA
!
      ARG2A = dCOS(c5)*dSIN(PLA90)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)
      BB    = dACOS(ARG2A)

      ARG2A = dSIN(c4)/dSIN(BB)
      ARG2A = dMAX1(ARG2A,-c1)
      ARG2A = dMIN1(ARG2A, c1)        
      ALPHA = dASIN(ARG2A)
!
!         Get GLOLIM - it is nesc. to test for the existence of solution
!
      ARGU2  = dCOS(c2)*dCOS(BB)/ &
                 (1.-dSIN(c4)*dSIN(BB)*dSIN(ALPHA))
      IF(dABS(ARGU2).GT.c1) THEN
      GLOLIM = 999.0
      ELSE
        GLOLIM = dACOS(ARGU2)/PI180
      END IF
!
!     Modify (if nesc.) the PLO solution
!
      IF((ABS(ARGU1).GT.GLOLIM.AND.GLAT.LE.RLAT).OR. &
        GLAT.GT.RLAT) THEN
            PLO1 = PI180*180.0 - PLO1
      END IF
!
!     The solution is symmetric so the direction must be if'ed
!
      IF(ARGU1.LT.0.0) THEN
         PLO1 = -PLO1
      END IF
!
!     Convert the radians to degrees
!
      PLA = PLA1/PI180        
      PLO = PLO1/PI180
!
!     To obtain a rotated value (ie so x-axis in pol.ste. points east)
!     add 270 to longitude
!
      PLO=AMOD(PLO+270.0,360.0)
!
 2000 CONTINUE      
      RETURN
      END                                  
!##################################################################
!
!      # Basic variables
!      #cname
!      #	   cvar
!      #	     cdimen
!      #		ccond1
!      #		  ccond2
!      #			 ctables-->
!      # Velocity, perturbation Exner function
!      UP	 : 1:3: 0:110000:vard:hist:anal:mpti:mpt3:mpt2
!      UC	 : 2:3: 0:110000:vard:hist:mpti:mpt3
!      VP	 : 3:3: 0:110000:vard:hist:anal:mpti:mpt3:mpt2
!      VC	 : 4:3: 0:110000:vard:hist:mpti:mpt3
!      WP	 : 5:3: 0:110000:vard:hist:anal:mpti:mpt3:mpt2
!      WC	 : 6:3: 0:110000:vard:hist:mpti:mpt3
!      PP	 : 7:3: 0:110000:vard:hist:anal:mpti:mpt3:mpt2
!      PC	 : 8:3: 0:110000:vard:hist:mpti:mpt3
!      # Theta_il, Hydrometeor mixing ratios
!      THP	 : 9:3: 0:110000:vard:sclp:hist:mpti:mpt3:mpt2
!      RTP	 :10:3: 1:110000:vard:sclp:hist:mpti:mpt3:mpt1
!      RCP	 :11:3: 2:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      RRP	 :12:3: 3:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      RPP	 :13:3: 4:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      RSP	 :14:3: 5:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      RAP	 :15:3: 6:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      RGP	 :16:3: 7:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      RHP	 :17:3: 8:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      # Hydrometeor concentrations
!      CCP	 :18:3: 9:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      CRP	 :19:3:10:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      CPP	 :20:3:11:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      CSP	 :21:3:12:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      CAP	 :22:3:13:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      CGP	 :23:3:14:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      CHP	 :24:3:15:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      # CCN and IN concentrations, pcp int energy, vapor, theta
!      CCCNP   :25:3:16:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      CCINP   :26:3:17:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      CDINP   :27:3:18:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      Q2	 :28:3: 3:110000:vard:hist:anal:mpti:mpt3
!      Q6	 :29:3: 7:110000:vard:hist:anal:mpti:mpt3
!      Q7	 :30:3: 8:110000:vard:hist:anal:mpti:mpt3
!      RV	 :31:3: 1:110000:vard:hist:anal:mpt3
!      THETA   :32:3: 0:110000:vard:hist:anal:mpti:mpt3
!      # TKE, cumulus parameterization, radiation
!      TKEP	 :33:3:19:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      #TKLP   :34:3: 0:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      THSRC   :35:3:20:110000:vard:hist:anal:mpti:mpt3
!      RTSRC   :36:3:21:110000:vard:hist:anal:mpti:mpt3
!      FTHRD   :37:3:22:110000:vard:hist:anal:mpti:mpt3
!      # Variable initialization
!      VARUP   :38:3:23:110000:vard:mpti
!      VARVP   :39:3:23:110000:vard:mpti
!      VARTP   :40:3:23:110000:vard:mpti
!      VARRP   :41:3:23:110000:vard:mpti
!      VARUF   :42:3:23:110000:vard:mpti
!      VARVF   :43:3:23:110000:vard:mpti
!      VARTF   :44:3:23:110000:vard:mpti
!      VARRF   :45:3:23:110000:vard:mpti
!      VARWTS  :46:3:23:110000:vard:mpti
!      # Reference state
!      PI0	 :47:3: 0:110000:vard:mpti
!      DN0	 :48:3: 0:110000:vard:mpti
!      TH0	 :49:3: 0:110000:vard:mpti
!      # Added scalars
!      SCLP	 :24:3: 0:110000:vard:sclp:hist:anal:mpti:mpt3:mpt1
!      # Tendencies
!      UT	 : 1:3: 0:100000:vari
!      VT	 : 2:3: 0:100000:vari
!      WT	 : 3:3: 0:100000:vari
!      PT	 : 4:3: 0:100000:vari
!      THT	 : 5:3: 0:100000:vari:sclt
!      RTT	 : 6:3: 1:100000:vari:sclt
!      RCT	 : 7:3: 2:100000:vari:sclt
!      RRT	 : 8:3: 3:100000:vari:sclt
!      RPT	 : 9:3: 4:100000:vari:sclt
!      RST	 :10:3: 5:100000:vari:sclt
!      RAT	 :11:3: 6:100000:vari:sclt
!      RGT	 :12:3: 7:100000:vari:sclt
!      RHT	 :13:3: 8:100000:vari:sclt
!      CCT	 :14:3: 9:100000:vari:sclt
!      CRT	 :15:3:10:100000:vari:sclt
!      CPT	 :16:3:11:100000:vari:sclt
!      CST	 :17:3:12:100000:vari:sclt
!      CAT	 :18:3:13:100000:vari:sclt
!      CGT	 :19:3:14:100000:vari:sclt
!      CHT	 :20:3:15:100000:vari:sclt
!      CCCNT   :21:3:16:100000:vari:sclt
!      CCINT   :22:3:17:100000:vari:sclt
!      CDINT   :23:3:18:100000:vari:sclt
!      TKET	 :24:3:19:100000:vari:sclt
!      #TKLT   :25:3: 0:100000:vari:sclt
!      # Scratch arrays
!      VT3DA   :26:3: 0:110000:vari
!      VT3DB   :27:3: 0:110000:vari
!      VT3DC   :28:3: 0:100000:vari
!      VT3DD   :29:3: 0:100000:vari
!      VT3DE   :30:3: 0:100000:vari
!      VT3DF   :31:3: 0:100000:vari
!      VT3DG   :32:3: 0:100000:vari
!      VT3DH   :33:3: 0:100000:vari
!      VT3DI   :34:3: 0:100000:vari
!      VT3DJ   :35:3: 0:100000:vari
!      VT3DK   :36:3: 0:100000:vari
!      VT3DL   :37:3: 0:100000:vari
!      VT3DM   :38:3: 0:100000:vari
!      VT3DN   :39:3: 0:100000:vari
!      VT3DO   :40:3: 0:100000:vari
!      VT3DP   :41:3: 0:100000:vari

!      # Tendencies of added scalars
!      SCLT	 :42:3:24:100000:vari:sclt

!      # Terrain-following coordinates and transformations
!      TOPT	 : 1:2: 0:110000:vard:hist:anal:mpti
!      TOPU	 : 2:2: 0:110000:vard:mpti
!      TOPV	 : 3:2: 0:110000:vard:mpti
!      TOPM	 : 4:2: 0:110000:vard:mpti
!      RTGT	 : 5:2: 0:110000:vard:mpti
!      RTGU	 : 6:2: 0:110000:vard:mpti
!      RTGV	 : 7:2: 0:110000:vard:mpti
!      RTGM	 : 8:2: 0:110000:vard:mpti
!      F13T	 : 9:2: 0:110000:vard:mpti
!      F13U	 :10:2: 0:110000:vard:mpti
!      F13V	 :11:2: 0:110000:vard:mpti
!      F13M	 :12:2: 0:110000:vard:mpti
!      F23T	 :13:2: 0:110000:vard:mpti
!      F23U	 :14:2: 0:110000:vard:mpti
!      F23V	 :15:2: 0:110000:vard:mpti
!      F23M	 :16:2: 0:110000:vard:mpti

!      # Horizontal coordinates and transformations
!      DXU	 :17:2: 0:110000:vard:mpti
!      DXV	 :18:2: 0:110000:vard:mpti
!      DXT	 :19:2: 0:110000:vard:mpti
!      DXM	 :20:2: 0:110000:vard:mpti
!      DYU	 :21:2: 0:110000:vard:mpti
!      DYV	 :22:2: 0:110000:vard:mpti
!      DYT	 :23:2: 0:110000:vard:mpti
!      DYM	 :24:2: 0:110000:vard:mpti
!      FMAPU   :25:2: 0:110000:vard:mpti
!      FMAPV   :26:2: 0:110000:vard:mpti
!      FMAPT   :27:2: 0:110000:vard:mpti
!      FMAPM   :28:2: 0:110000:vard:mpti
!      FMAPUI  :29:2: 0:110000:vard:mpti
!      FMAPVI  :30:2: 0:110000:vard:mpti
!      FMAPTI  :31:2: 0:110000:vard:mpti
!      FMAPMI  :32:2: 0:110000:vard:mpti
!      GLAT	 :33:2: 0:110000:vard:mpti
!      GLON	 :34:2: 0:110000:vard:mpti

!      # Surface fluxes
!      UW	 :35:2: 0:110000:vard:anal:mpt3
!      VW	 :36:2: 0:110000:vard:anal:mpt3
!      WFZ	 :37:2: 0:110000:vard:anal:mpt3
!      TFZ	 :38:2: 0:110000:vard:anal:mpt3
!      QFZ	 :39:2: 0:110000:vard:anal:mpt3

!      # Surface precipitation
!      ACCPR   :40:2: 0:110000:vard:hist:anal:mpti:mpt3
!      ACCPP   :41:2: 0:110000:vard:hist:anal:mpti:mpt3
!      ACCPS   :42:2: 0:110000:vard:hist:anal:mpti:mpt3
!      ACCPA   :43:2: 0:110000:vard:hist:anal:mpti:mpt3
!      ACCPG   :44:2: 0:110000:vard:hist:anal:mpti:mpt3
!      ACCPH   :45:2: 0:110000:vard:hist:anal:mpti:mpt3
!      PCPRR   :46:2: 0:110000:vard:anal:mpt3
!      PCPRP   :47:2: 0:110000:vard:anal:mpt3
!      PCPRS   :48:2: 0:110000:vard:anal:mpt3
!      PCPRA   :49:2: 0:110000:vard:anal:mpt3
!      PCPRG   :50:2: 0:110000:vard:anal:mpt3
!      PCPRH   :51:2: 0:110000:vard:anal:mpt3
!      PCPG    :52:2: 0:110000:vard:hist:anal:mpti:mpt3
!      QPCPG   :53:2: 0:110000:vard:hist:anal:mpti:mpt3
!      DPCPG   :54:2: 0:110000:vard:hist:anal:mpti:mpt3
!      ACONPR  :55:2: 0:110000:vard:hist:anal:mpti:mpt3
!      CONPRR  :56:2: 0:110000:vard:hist:anal:mpt3

!      # Surface radiation
!      RSHORT  :57:2: 0:110000:vard:hist:anal:mpti:mpt3
!      RLONG   :58:2: 0:110000:vard:hist:anal:mpti:mpt3
!      RLONGUP :59:2: 0:110000:vard:hist:anal:mpti:mpt3
!      ALBEDT  :60:2: 0:110000:vard:hist:anal:mpti:mpt3

!      # Variable initialization, Coriolis
!      VARP	 :61:2: 0:110000:vard:mpti
!      VARP2   :62:2: 0:110000:vard:mpti
!      FCORU   :63:2: 0:110000:vard:mpti
!      FCORV   :64:2: 0:110000:vard:mpti

!      # Scratch arrays
!      VT2DA   :65:2: 0:110000:vard
!      VT2DB   :66:2: 0:100000:vard
!      VT2DC   :67:2: 0:100000:vard
!      VT2DD   :68:2: 0:100000:vard
!      VT2DE   :69:2: 0:100000:vard
!      VT2DF   :70:2: 0:100000:vard

!      # LEAF-2 4-D arrays
!      TGP	 : 1:s: 0:110000:vard:hist:anal:mpti:mpt3
!      WGP	 : 2:s: 0:110000:vard:hist:anal:mpti:mpt3
!      SCHAR     : 3:s: 0:110000:vard:hist:anal:mpti:mpt3
!      GSF	 : 4:s: 0:110000:vard:hist:anal:mpti:mpt3

!      # Tables:
!      #--------

!      # vard - grid dependent variable
!      # vari - grid independent variable
!      # hist - write to history file
!      # anal - write to analysis file
!      # sclp - scalar table - past time level
!      # sclt - scalar table - tendency


!      # Parallel tables:
!      #-----------------

!      #  mpti - initialization, full sub-domain master to node
!      #  mpt1 - long timestep, subdomain boundaries node to node
!      #  mpt2 - small timestep, subdomain boundaries node to node
!      #  mpt3 - full sub-domain node to master for output
!      # description
!      #   1      2 3  4   5     6 +
!      #-----------------------------------------------------
!      #  RTP   : 8:3: 3:110000:vard:sclp:hist:mpti:mpt3:mpt1
!      # 1. variable name
!      # 2. variable number (from tables at start)
!      # 3. dimensionality ( 3-> 3d; 2-> 2d; s->soil)
!      # 4. conditional number
!      # 5. option flags :110000:
!      #		    1- node gets this variable? (0-no, 1-yes)
!      #		     2- master gets this variable? (0-no, 1-yes)
!      #		      3- not used
!      #		       4- not used
!      #			5- not used
!      #			 6- not used
!      # 6.+ list of tables
!#########################################################################

SUBROUTINE HTINT(NZZ1,VCTRA,ELEVA,NZZ2,VCTRB,ELEVB)
DIMENSION VCTRA(*),VCTRB(*),ELEVA(*),ELEVB(*)

L=1
DO 20 K=1,NZZ2
30 CONTINUE
IF(ELEVB(K).LT.ELEVA(1))GO TO 35
IF(ELEVB(K).GE.ELEVA(L).AND.ELEVB(K).LE.ELEVA(L+1))GO TO 35
IF(ELEVB(K).GT.ELEVA(NZZ1))GO TO 36
L=L+1
IF(L.EQ.NZZ1) then
  print *,'htint:nzz1',nzz1
  do kk=1,L
    print*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
  enddo
  stop 'htint'
endif
GO TO 30
35 CONTINUE
WT=(ELEVB(K)-ELEVA(L))/(ELEVA(L+1)-ELEVA(L))
VCTRB(K)=VCTRA(L)+(VCTRA(L+1)-VCTRA(L))*WT
GO TO 20
36 CONTINUE
WT=(ELEVB(K)-ELEVA(NZZ1))/(ELEVA(NZZ1-1)-ELEVA(NZZ1))
VCTRB(K)=VCTRA(NZZ1)+(VCTRA(NZZ1-1)-VCTRA(NZZ1))*WT
20 CONTINUE

return 
end

!SUBROUTINE AE0(NPTS,A,B)
!DIMENSION A(NPTS)
!DO I=1,NPTS
!  A(I)=B
!ENDDO
!RETURN
!END


FUNCTION RS(P,T)

ES=610.78*EXP(17.269*(T-273.16)/(T-35.86))
RS=.622*ES/(P-ES)

RETURN
END

function td(p,rs)

implicit none
real rr,rs,es,esln,p,td

rr=rs+1e-8
es=p*rr/(.622+rr)
esln=log(es)
td=(35.86*esln-4947.2325)/(esln-23.6837)

return
end
