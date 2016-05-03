cc This is a subroutine which computes the nearest neighbor list based on user input
cc of the radius cut-off value 'rcut' and its squared value 'rcut2'.
cc The routines which call this routine have not been included, so that it may be
cc inserted into a variety of different codes with ease.
cc
cc Here is a list of the terminology and definitions:
c perlen(i) is the simulation cell length in i dimeion
c numnbr(i) is the number of neighbors for atom i
c natoms is the number of atoms, read in from dump file
c neimax is the maximum number of neighbors, predetermined
c disnbr(a,b,c) is the distance between atom b and c in the a dimension
c Rf(i,j) is the reference position of atom j, for the i dimension
c reflb(i) is the reference lower bound of i dimension
c refub(i) is the reference upper bound of i dimension
c reflen(i) is the reference length of i dimension
c inbr(i,j) is the atom number of neighbor i of atom j
c
cc the first 3 places in the disnbr(a,j,i) array stores the curent distance between
cc atom i and neighbor j, while the last 3 places keep track of the original
cc distance between atom i and neighbor j.
cc
cc The calculation of interatomic distances in this code assumes periodic boundary
cc conditions. If for some reason, this is not valid, the user will need to
cc modify the lines of the code that use 'nint' which deals with nearest image calculations.
cc
c*************************************************************************
      subroutine setnbr
	  use data
      implicit real*4(a-h,o-z)
      implicit integer(i-n)
      dimension per12(3),dis(6)
      dimension bndlnk(2,3),dlnk(3),nglnk(3),lnklst(0:10*natoms)
      dimension nltmp(10*neimax),isper(3)
      data icall/0/ ! which means icall = 0
      save bndlnk,dlnk,nglnk,isper ! sets constant after procedure call
      real*8 dist1,dist2
      integer atoms
c there are two updates performed:
c one is for each call to setnbr,
c for this one we update r1; this makes the calculation of disnbr easy
c the other is for when some atom moves enough, then we update r0
c on the first time through, we calculate everything
c
c Keep a list of atom positions since last update of neighbor list.
c If any atom has moved more than drcut/2., then update neighbor list.
c Otherwise, keep same list.

c
c check displacements
      per12(1) = 0.5*perlen(1)
      per12(2) = 0.5*perlen(2)
      per12(3) = 0.5*perlen(3)
c Update numbnr, inbr, and disnbr.
      nupdates = nupdates + 1
c update neighbor list
      rbuf = rcut
      rbuf2 = rbuf**2
      rcut2 = rcut**2
      do i=1,natoms
        numnbr(i) = 0
      end do
      do i=1,natoms
         do j=1,neimax
            disnbr(1,j,i) = 0.0
            disnbr(2,j,i) = 0.0
            disnbr(3,j,i) = 0.0
            disnbr(4,j,i) = 0.0
            disnbr(5,j,i) = 0.0
            disnbr(6,j,i) = 0.0
         enddo
      enddo
c Storage of neighbor indices with determination of neighbor list
c via linked lists.
c determine the upper and lower bounds of the atoms coordinates
      do 4110 kc = 1,3
        bndlnk(1,kc) = 1.e6
        bndlnk(2,kc) = -1.e6
c         bndlnk(1,kc) = perlb(kc)
c	   bndlnk(2,kc) = perub(kc)
4110  continue
      do 4120 j = 1,natoms
        do 4125 kc = 1,3
          rtmp = Rf(kc,j)
          if (rtmp.gt.refub(kc)) then
            rtmp = rtmp - reflen(kc)*
     &             float(1 + int((rtmp-refub(kc))/reflen(kc)))
          else if (rtmp.lt.reflb(kc)) then
            rtmp = rtmp + reflen(kc)*
     &             float(1 + int((reflb(kc)-rtmp)/reflen(kc)))
          endif
          bndlnk(1,kc) = min(bndlnk(1,kc),rtmp)
          bndlnk(2,kc) = max(bndlnk(2,kc),rtmp)
4125    continue
4120  continue
      do 4127 kc = 1,3
        isper(kc) = 1
        if (abs(bndlnk(1,kc)-reflb(kc)).lt.5.) then
          bndlnk(1,kc) = reflb(kc)
        else
          bndlnk(1,kc) = bndlnk(1,kc) - 1.0
          isper(kc) = 0
        endif
        if (abs(refub(kc) - bndlnk(2,kc)).lt.5.) then
          bndlnk(2,kc) = refub(kc)
        else
          bndlnk(2,kc) = bndlnk(2,kc) + 1.0
          isper(kc) = 0
        endif
4127  continue
      voltmp = (bndlnk(2,1)-bndlnk(1,1))
     &        *(bndlnk(2,2)-bndlnk(1,2))
     &        *(bndlnk(2,3)-bndlnk(1,3))
      write(6,*) 'natoms=',natoms
      voltmp = voltmp/(10.*natoms)
      dltmp = voltmp**(1./3.)
      write(6,*) 'dltmp=',dltmp
      do 4130 kc = 1,3
c nglnk determines the box count in the kc dimension
        nglnk(kc) = int((bndlnk(2,kc)-bndlnk(1,kc))/dltmp)
        nglnk(kc) = max(3,nglnk(kc))
c dlnk determines the box size in the kc dimension
        dlnk(kc) = (bndlnk(2,kc)-bndlnk(1,kc))/nglnk(kc)
4130  continue
c clear the link list array
c ngtot gives the total box count in the simulation domain
      ngtot = nglnk(1)*nglnk(2)*nglnk(3)
      if (ngtot.gt.10*natoms) stop 'gneigh 4'
c clear each boxcount
      do 4140 kc = 0,ngtot
        lnklst(kc) = 0
4140  continue

c place the atoms in the linked list
      do 4200 j = 1,natoms
        ix = intmod(int((Rf(1,j)-bndlnk(1,1))/dlnk(1)),
     &       nglnk(1))
        iy = intmod(int((Rf(2,j)-bndlnk(1,2))/dlnk(2)),
     &       nglnk(2))
        iz = intmod(int((Rf(3,j)-bndlnk(1,3))/dlnk(3)),
     &       nglnk(3))
        il = ix*nglnk(2)*nglnk(3) + iy*nglnk(3) + iz
c lnklst(il) will be the box where atom j is located
        if (lnklst(il).eq.0) then
          lnklst(il) = j
        else
          print *,"two atoms in one cell of linked list"
          print *,il,lnklst(il),j,ix,iy,iz
          print *,(Rf(kc,lnklst(il)),kc=1,3)
          print *,(Rf(kc,j),kc=1,3)
c          stop
        endif
4200  continue
c     write(6,*) 'isper:',isper
c     write(6,9096) 'bndlnk:',bndlnk
9096  format(a,2x,6f9.3)
c     write(6,*) 'nglnk:',nglnk
c     write(6,*) 'dlnk:',dlnk
c
c determine new neighbor list while getting the neighbors
c
      do 4000 i=1,natoms
c       nnindx(i) = nnindx(i-1)
        nnltmp = 0
        ix0 = intmod(int((Rf(1,i)-bndlnk(1,1))/dlnk(1)),
     &        nglnk(1))
        iy0 = intmod(int((Rf(2,i)-bndlnk(1,2))/dlnk(2)),
     &        nglnk(2))
        iz0 = intmod(int((Rf(3,i)-bndlnk(1,3))/dlnk(3)),
     &        nglnk(3))
        drtmp = sqrt(rbuf2)
c ndx is the box count version of rbuf2**.5
        ndx = 1 + int(drtmp/dlnk(1))
        if (isper(1).eq.0) then
          ndx1 = max(0,ix0-ndx)
          ndx2 = min(nglnk(1)-1,ix0+ndx)
        else
          ndx1 = ix0 - ndx
          ndx2 = ix0 + ndx
cc        ndx2 = ndx1 + min(2*ndx,nglnk(1)-1)
        endif
        ndy = 1 + int(drtmp/dlnk(2))
        if (isper(2).eq.0) then
          ndy1 = max(0,iy0-ndy)
          ndy2 = min(nglnk(2)-1,iy0+ndy)
        else
          ndy1 = iy0 - ndy
          ndy2 = iy0 + ndy
cc        ndy2 = ndy1 + min(2*ndy,nglnk(2)-1)
        endif
        ndz = 1 + int(drtmp/dlnk(3))
        if (isper(3).eq.0) then
          ndz1 = max(0,iz0-ndz)
          ndz2 = min(nglnk(3)-1,iz0+ndz)
        else
          ndz1 = iz0 - ndz
          ndz2 = iz0 + ndz
cc        ndz2 = ndz1 + min(2*ndz,nglnk(3)-1)
        endif
c
c Check if atom is close to the periodic boundaries. If not, do the first
c set of loops here to save CPU time.  Only do the second set of loops
c for atoms near the edges.
c (Idea from John Hamilton, 1/30/97, who tested it for 2D.)
c
        if(ndx1.le.0) go to 4505
        if(ndy1.le.0) go to 4505
        if(ndz1.le.0) go to 4505
        if(ndx2.ge.nglnk(1)) go to 4505
        if(ndy2.ge.nglnk(2)) go to 4505
        if(ndz2.ge.nglnk(3)) go to 4505
c
        do ix = ndx1,ndx2
          do iy = ndy1,ndy2
            do iz = ndz1,ndz2
              il = ix*nglnk(2)*nglnk(3) + iy*nglnk(3) + iz
              j = lnklst(il)
              if (j.ne.0.and.j.lt.i) then
                nnltmp = nnltmp + 1
                if (nnltmp.gt.10*neimax) then
                  write(6,*)' i, j, nnltmp = ',i,j,nnltmp
                  stop 'nnltmp 1'
                end if
                nltmp(nnltmp) = j
              endif
            end do
          end do
        end do
        go to 4540
c
4505    continue
c Reset loop ranges and use intmod.
        if (isper(1).ne.0) ndx2 = ndx1 + min(2*ndx,nglnk(1)-1)
        if (isper(2).ne.0) ndy2 = ndy1 + min(2*ndy,nglnk(2)-1)
        if (isper(3).ne.0) ndz2 = ndz1 + min(2*ndz,nglnk(3)-1)
        do 4510 ixt = ndx1,ndx2
          ix = intmod(ixt,nglnk(1))
          do 4520 iyt = ndy1,ndy2
            iy = intmod(iyt,nglnk(2))
            do 4530 izt = ndz1,ndz2
              iz = intmod(izt,nglnk(3))
              il = ix*nglnk(2)*nglnk(3) + iy*nglnk(3) + iz
              j = lnklst(il)
              if (j.ne.0.and.j.lt.i) then
                nnltmp = nnltmp + 1
                if (nnltmp.gt.10*neimax) then
                  write(6,*)' i, j, nnltmp = ',i,j,nnltmp
                  stop 'nnltmp 2'
                end if
                nltmp(nnltmp) = j
              endif
4530        continue
4520      continue
4510    continue
c
4540    continue
c
        call isort(nnltmp,nltmp)
c nnltmp is the nearest neighbor count for atom i
        do 4600 jj = 1,nnltmp
c nltmp gives the tag of atom is nearest neighbor j
          j = nltmp(jj)
c compute the square of the distance to the closest periodic image
          dis(4) = Rf(1,i) - Rf(1,j)
          dis(4) = dis(4) - reflen(1)*nint(dis(4)/reflen(1))
          dis2 = dis(4)**2
          if (dis2.gt.rbuf2) goto 4600
          dis(5) = Rf(2,i) - Rf(2,j)
          dis(5) = dis(5) - reflen(2)*nint(dis(5)/reflen(2))
          dis2 = dis2 + dis(5)**2
          if (dis2.gt.rbuf2) goto 4600
          dis(6) = Rf(3,i) - Rf(3,j)
          dis(6) = dis(6) - reflen(3)*nint(dis(6)/reflen(3))
          dis2 = dis2 + dis(6)**2
c
c determine if these particles are within the storage distance
          if (dis2.gt.rbuf2) goto 4600
c
c store the index of the particle
c         nnindx(i) = nnindx(i) + 1
c         nnlst(nnindx(i)) = j
c
c determine which pairs are separated by less than rcut
c and store the needed information about these pairs
c
c if nearest periodic image is out of range, then all images will be
c
          if (dis2.gt.rcut2) go to 4600
c
          numnbr(i) = numnbr(i) + 1
          numnbr(j) = numnbr(j) + 1
          ni = numnbr(i)
          nj = numnbr(j)
          inbr(ni,i) = j
          inbr(nj,j) = i
c Needed to switch the signs here relative to the original version of
c kroner3.f, since i and j are defined oppositely in gneigh (dynamo code)
c where the linked neighbor list code was pulled from.
          disnbr(4,ni,i) = -dis(4)
          disnbr(5,ni,i) = -dis(5)
          disnbr(6,ni,i) = -dis(6)
          disnbr(4,nj,j) = dis(4)
          disnbr(5,nj,j) = dis(5)
          disnbr(6,nj,j) = dis(6)
c
          dis(1) = r(1,i) - r(1,j)
          dis(1) = dis(1) - perlen(1)*nint(dis(1)/perlen(1))
          dis(2) = r(2,i) - r(2,j)
          dis(2) = dis(2) - perlen(2)*nint(dis(2)/perlen(2))
          dis(3) = r(3,i) - r(3,j)
          dis(3) = dis(3) - perlen(3)*nint(dis(3)/perlen(3))

          disnbr(1,ni,i) = -dis(1)
          disnbr(2,ni,i) = -dis(2)
          disnbr(3,ni,i) = -dis(3)
          disnbr(1,nj,j) = dis(1)
          disnbr(2,nj,j) = dis(2)
          disnbr(3,nj,j) = dis(3)

c
cc Assume that perlen(3) > rbuf so no images in z-direction.
c
4600    continue
4000  continue
9099      format(6x,i6,3f9.3,i6)
      nsum=0
      maxnbr=0
      do 1200 i=1,natoms
c       write(6,*)'i, numnbr(i) = ',i,numnbr(i)
        nsum=nsum+numnbr(i)
        maxnbr=max(maxnbr,numnbr(i))
        if(numnbr(i).gt.neimax)then
          write(6,9090)i,numnbr(i)
9090      format('  atom ',i5,' has ',i5,' neighbors',/,
     &           ' exceeds dimension ')
          stop
        endif
1200  continue
      ave=real(nsum)/real(natoms)
c     write(6,*)'maximum number of neighbors = ',maxnbr
c     write(6,*)'average number of neighbors = ',ave
c
      return
      end
c
c************************************************************************
      subroutine isort(n,ra)
        use data

      implicit real*4(a-h,o-z)
      implicit integer(i-n)
      integer ra(n),rra
c smf: if 0 or 1 element in list just stop here because otherwise the
c code blows up and there is nothing to do anyway
      if (n.lt.2) return
      l=n/2+1
      ir=n
   10 continue
      if(l.gt.1)then
         l=l-1
         rra=ra(l)
      else
         rra=ra(ir)
         ra(ir)=ra(1)
         ir=ir-1
         if(ir.eq.1)then
            ra(1)=rra
            return
         endif
      endif
      i=l
      j=l+l
   20 if(j.le.ir)then
         if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
         endif
         if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
         else
            j=ir+1
         endif
         go to 20
      endif
      ra(i)=rra
      go to 10
      end
c ***************************************************************************
      subroutine asort(n,ra,indxra)
        use data

      implicit real*4(a-h,o-z)
      implicit integer(i-n)
      dimension ra(n)
      integer indxra(2,n)
c smf: if 0 or 1 element in list just stop here because otherwise the
c code blows up and there is nothing to do anyway
      if (n.lt.2) return
      l=n/2+1
      ir=n
   10 continue
      if(l.gt.1)then
         l=l-1
         rra=ra(l)
         ii1=indxra(1,l)
         ii2=indxra(2,l)
      else
         rra=ra(ir)
         ii1=indxra(1,ir)
         ii2=indxra(2,ir)
         ra(ir)=ra(1)
         indxra(1,ir)=indxra(1,1)
         indxra(2,ir)=indxra(2,1)
         ir=ir-1
         if(ir.eq.1)then
            ra(1)=rra
            indxra(1,1)=ii1
            indxra(2,1)=ii2
            return
         endif
      endif
      i=l
      j=l+l
   20 if(j.le.ir)then
         if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
         endif
         if(rra.lt.ra(j))then
            ra(i)=ra(j)
            indxra(1,i)=indxra(1,j)
            indxra(2,i)=indxra(2,j)
            i=j
            j=j+j
         else
            j=ir+1
         endif
         go to 20
      endif
      ra(i)=rra
      indxra(1,i)=ii1
      indxra(2,i)=ii2
      go to 10
      end
c
c************************************************************************
      function intmod(i,n)
      implicit real*4 (a-h,o-z)
      if (i.ge.0) then
         intmod = mod(i,n)
      else
         intmod = n - 1 + mod(i+1,n)
      endif
      return
      end
c
c************************************************************************
      function cmpsum(i)

	use data
      implicit real*4 (a-h,o-z)
c      parameter (natmax=857601,neimax=30)
c      common /lattice/ r(3,natmax),Rf(3,natmax),itype(natmax),outflag,
c     &                 reflen(3),refub(3),reflb(3),csym(natmax),
c     &                 v(3,natmax),
c     &                 eng(natmax),cd(natmax),s(6,natmax),
c     &                 W(3,3,natmax),Phi(3,natmax),gradPhi(3,3,natmax),
c     &                 Vorticity(3,3,natmax),Strain(3,3,natmax)
c     &                 ,svalue(natmax)
c      common /neigh/ numnbr(natmax),inbr(neimax,natmax),
c     &               disnbr(6,neimax,natmax)
c	common /neigh_new/ numnbr_new(natmax),inbr_new(neimax,natmax),
c     &               disnbr_new(3,neimax,natmax)
c      common /vecs/ alat,vecsum(3,natmax)
      dimension sumtmp(neimax*(neimax-1)/2)
      integer indxst(2,neimax*(neimax-1)/2)
      dimension term(3)
c
      btol = 0.5*(alat/sqrt(6.0))
      npairs = 0
      icnt = 0
      do k = 1,3
         vecsum(k,i) = 0.0
      enddo
      do i1 = 1,numnbr(i)-1
         do i2 = i1+1,numnbr(i)
            npairs = npairs + 1
            indxst(1,npairs) = i1
            indxst(2,npairs) = i2
            sumtmp(npairs) = 0.0
            do k = 1,3
               sumtmp(npairs) = sumtmp(npairs)
     $              + (disnbr(k,i1,i)+disnbr(k,i2,i))**2
            enddo
         enddo
         termmag = 0.0
         do k = 1,3
            term(k)=disnbr(k,i1,i)-disnbr(k+3,i1,i)
            vecsum(k,i) = vecsum(k,i) - term(k)
            termmag = termmag + term(k)**2
         enddo
         if (sqrt(termmag).gt.btol) icnt = icnt + 1
      enddo
c
      if (icnt.gt.0) then
         vs = 0.0
         do k = 1,3
            vecsum(k,i) = vecsum(k,i) / icnt
            vs = vs + vecsum(k,i)**2
         enddo
c        if (vs.gt.0.1) print*,i,sqrt(vs),icnt
      endif
c
      call asort(npairs,sumtmp,indxst)
      sum = 0.0
      do i1 = 1,6
         sum = sum + sumtmp(i1)
      enddo
      cmpsum = sum
      return
      end
c*************************************************************************
      subroutine scale(vec,sca)
        use data

      implicit real*4(a-h,o-z)
      implicit integer(i-n)
      dimension vec(1)
      do 10 i=1,3
10      vec(i) = vec(i)*sca
      return
      end
c*************************************************************************
      function cvmgp(x1,x2,x3)
      implicit real*4 (a-h,o-z)
c
      if(x3.gt.0.0)then
       cvmgp=x1
      else
       cvmgp=x2
      endif
c     write(66,'(a,4f15.5)')'x1,x2,x3,cvmgp =',x1,x2,x3,cvmgp
      return
      end
c
c**********************************************************************
