C -----------------------------------------------------------C
C Subroutine for Abaqus/Explicit for isotropic elasticity    C
C and isotropic plasticity with pressure and lode            C
C dependence, according to plasticity model by BAI, Y. and   C
C WIERZBICKI, T. (published in "A new model of metal         C
C plasticity and fracture with pressure and Lode             C
C dependence", International Journal of Plasticity 24 (2008) C
C pages 1071-1096).                                          C
C															 C
C - Damage initiation according to 3D fracture locus         C
C   (equivalent plastic strain to failure initiation being   C
C   function of both hydrostatic pressure and lode angle)    C
C - Damage evolution based on energy dissipated during       C
C   damage process                                           C
C															 C
C - Temperature softening									 C
C - Strain rate hardening									 C
C                                                            C
C Originally based on the von Mises' J2 plasticity theory.   C
C Elastic predictor, radial corrector algorithm used. With   C
C explicit forward Euler scheme for integration of flow rule.C
C                                                            C
C Not suitable for solving problem settings including plane  C
C stresses. Not suitable for 2D models.                      C
C -----------------------------------------------------------C
C Solution dependent variables (SDVs):						 C
C															 C
C SDV1: 	Equivalent plastic strain						 C
C SDV2: 	Damage variable									 C
C SDV3: 	Yield stress at damage initiation				 C
C SDV4: 	Flag (=0 element not damaged,					 C
C			      =1 element experienced brittle damage,	 C
C			      =2 element experienced ductile damage)	 C
C SDV5-10:	Total strain tensor								 C
C SDV11-16:	Plastic strain tensor							 C
C SDV17:	Equivalent plastic strain increment				 C
C SDV18:	Temperature softening correction function		 C
C SDV19:	Strain rate hardening correction function		 C
C SDV20:	Equivalent plastic strain rate					 C
C SDV21:	Equivalent plastic strain rate (computed over	 C
C           a user-defined no. of increments)				 C
C SDV22:	Increment counter (for eq. plas. str. rate       C
C           computation)									 C
C SDV23:	Sum of equivalent plastic strain increments      C
C           before the computation of eq. plas. str. rate    C
C SDV24:	Sum of time increments before computation of eq. C
C           plas. str. rate									 C
C SDV25:	vacant											 C
C SDV26:	Maximum principal stress						 C
C SDV27:	Equivalent plastic strain (averaged over a		 C
C			user-defined no. of increments)					 C
C SDV28:	Element deletion flag							 C
C -----------------------------------------------------------C
C Contact info:                                              C
C															 C
C Mohamed Sharaf                                             C
C RWTH Aachen                                                C
C Institut fuer Eisenhuettenkunde                            C
C Intzestrasse 1                                             C
C 52072 Aachen                                               C
C Germany                                                    C
C mohamed.sharaf@iehk.rwth-aachen.de                         C
C                                                            C
C Aachen, 6 March 2012                                       C
C -----------------------------------------------------------C
        subroutine vumat(
C Read only
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     1  stepTime, totalTime, dt, cmname, coordMp, charLength,
     1  props, density, strainInc, relSpinInc,
     1  tempOld, stretchOld, defgradOld, fieldOld,
     1  stressOld, stateOld, enerInternOld, enerInelasOld,
     1  tempNew, stretchNew, defgradNew, fieldNew,
C Write only
     1  stressNew, stateNew, enerInternNew, enerInelasNew)
C
        include 'vaba_param.inc'
C
        dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     1  relSpinInc(*), tempOld(nblock),
     1  stretchOld(*), defgradOld(*),
     1  fieldOld(*), stressOld(nblock,ndir+nshr),
     1  stateOld(nblock,nstatev), enerInternOld(nblock),
     1  enerInelasOld(nblock), tempNew(nblock),
     1  stretchNew(*), defgradNew(*), fieldNew(*),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     1  enerInternNew(nblock), enerInelasNew(nblock)
C
        character*80 cmname
C
C Defining numerical constants
        parameter( zero = 0., one = 1., two = 2., three = 3.,
     1  third = one/three, half = .5, twoThirds = two/three,
     1  threeHalves = 1.5, thousand = 1000., tol=1.D-6 )
C
C
C
C Defining material properties constants
        e      = props(1)
        xnu    = props(2)
        ceta   = props(3)
        eta0   = props(4)
        cthetas= props(5)
        cthetat= props(6)
        cthetac= props(7)
        om     = props(8)
C       Ductile damage initiation locus function parameters
        d1     = props(9)
        d2     = props(10)
        d3     = props(11)
        d4     = props(12)
        d5     = props(13)
        d6     = props(14)
C       Fracture energy dissipation (energy dissipated per element unit area after
C       ductile damage dissipation)
        gf     = props(15)
C       Maximum value for ductile damage variable
        Dcrit   = props(16)
C       5 parameters for temperature softening correction function
        cT1    = props(17)
        cT2    = props(18)
		cT3    = props(19)
        eta2   = props(20)
        cp     = props(21)
        t0     = props(22)
C       2 parameters for strain rate hardening correction function
        cE1    = props(23)
        cE2    = props(24)
		cE3    = props(25)
C       Number of increments to update plastic strain rate
        strrInc= props(26)
C       Brittle damage initiation stress value
        sigdmg = props(27)
C
        nvalue = (nprops-32)/2
C
C Checking for valid entries
        if (om.lt.zero) then
          write(6,5)
 5        format(//,30X,'***ERROR - m MUST BE A NON-NEGATIVE INTEGER')
        endif
C
C Defining constants
        twomu  = e / ( one + xnu )
        thremu = threeHalves * twomu
        sixmu  = three * twomu
        alamda = twomu * ( e - twomu ) / ( sixmu - two * e )
        akappa = e * twomu * half / (three * (threeHalves * twomu - e ))
        term   = one / ( twomu * ( one + hard/thremu ) )
        con1   = sqrt( twoThirds )
        pi     = 4. * atan(one)
        con2   = (sqrt(three)/two) / (one - sqrt(three)/two)
C
C Computation per material point starts here
      do 100 i = 1,nblock
C
C If brittle damage has previously occured, stress is zeroed
        if (stateOld(i,4).gt.half .and. stateOld(i,4).lt.(one+half)) then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
          stateNew(i,4) = one
          goto 100
        endif
C Updating total strain (states 5-10)
        stateNew(i,5)  = stateOld(i,5)  + strainInc(i,1)
        stateNew(i,6)  = stateOld(i,6)  + strainInc(i,2)
        stateNew(i,7)  = stateOld(i,7)  + strainInc(i,3)
        stateNew(i,8)  = stateOld(i,8)  + strainInc(i,4)
        stateNew(i,9)  = stateOld(i,9)  + strainInc(i,5)
        stateNew(i,10) = stateOld(i,10) + strainInc(i,6)
C
C       Trace of total strain tensor
        epsilontrace = stateNew(i,5) + stateNew(i,6) + stateNew(i,7)
C
C       Trace of strain increment tensor
        trace = strainInc(i,1) + strainInc(i,2) + strainInc(i,3)
C
C Calculating softening correction term due to temperature rise
        if (tempOld(i).le.tol .and. tempOld(i).ge.-tol) then
		  facT = cT1*exp(zero-cT2*t0)+cT3
          stateNew(i,18)= facT
        else
          facT = cT1*exp(zero-cT2*tempOld(i))+cT3
          stateNew(i,18)= facT
        endif   
C
C Calculating hardening correction term due to straining rate
		facStrrate = cE1*log(stateOld(i,21)) + cE2 + cE3*stateOld(i,21)
C
C Preventing negative values for straining rate
        if (facStrrate.lt. one) facStrrate = one
        stateNew(i,19)= facStrrate
C
C Restoring previous softening correction term due to ductile damage
        facDctlDmgDev = one - stateOld(i,2)
        facDctlDmgVol = facDctlDmgDev
        if (stressOld(i,1)+stressOld(i,2)+stressOld(i,3).lt.zero) then
          facDctlDmgVol = one
        endif
C
C Trial stress
        sig1 = stressOld(i,1) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,1)
        sig2 = stressOld(i,2) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,2)
        sig3 = stressOld(i,3) + facDctlDmgVol*alamda*trace
     1                        + facDctlDmgDev*twomu*strainInc(i,3)
        sig4 = stressOld(i,4) + facDctlDmgDev*twomu*strainInc(i,4)
        sig5 = stressOld(i,5) + facDctlDmgDev*twomu*strainInc(i,5)
        sig6 = stressOld(i,6) + facDctlDmgDev*twomu*strainInc(i,6)
C
C Following equation numbering of publication mentioned in code title text
C Equ. (4) - Calculating the deviatoric part of trial stress
        sigmean = third * ( sig1 + sig2 + sig3 )
        ds1 = sig1 - sigmean
        ds2 = sig2 - sigmean
        ds3 = sig3 - sigmean
C
C Calculating the magnitude of the deviatoric trial stress tensor
        dsmag = sqrt( ds1**2 + ds2**2 + ds3**2 + two*sig4**2 + two*sig5**2
     1  + two*sig6**2 )
C
C Preventing a divide by zero when dsmag is zero. When dsmag is zero, computation is still in the elastic zone
        if (dsmag.lt.tol .and. dsmag.ge.zero) dsmag = tol
        if (dsmag.gt.(zero-tol) .and. dsmag.le.zero) dsmag = zero-tol
C
C Following equation numbering of publication mentioned in code title text
C Eq. (1) - Calculating the 1st invariant of the stress tensor
        p   = zero - sigmean
C
C Eq. (2) - Calculating the 2nd invariant of the stress tensor
        q   = dsmag / con1
C
C Eq. (3) - Calculating the 3rd invariant of the stress tensor
        r   = ( three/twoThirds * (ds1*(ds1**2 + sig4**2 + sig6**2) + 2.
     1  *sig4*(ds1*sig4 + ds2*sig4 + sig6*sig5) + 2.*sig6*(ds1*sig6 + sig4
     1  *sig5 + ds3*sig6) + ds2*(sig4**2 + ds2**2 + sig5**2) + 2.*sig5
     1  *(sig4*sig6 + ds2*sig5 + ds3*sig5) + ds3*(sig6**2 + sig5**2
     1  + ds3**2)) )**third
C
C Eq. (5) - Calculating the dimensionless hydrostatic pressure eta
        eta = sigmean / q
C
C Eq. (6) - Calculating the Lode angle theta
        cosine=(r/q)**3
C       Assuring that -1<cosine<1
        if (cosine.gt.one) cosine=one
        if (cosine.lt.(zero-one)) cosine=zero-one
        theta = third*acos(cosine)
C
C Eq. (7) - Calculating the normalized Lode angle thetabar
        thetabar = one - 6.*theta/pi
C
C Eq. (12) - Calculating the Lode angle dependant parameter gamma
        gamma = con2 * (one/cos(theta-pi/6.) - one)
C
C Eq. (13) - Determining cthetaax, either tension case or compression case
        if( thetabar .ge. zero ) then
          cthetaax = cthetat
        else
          cthetaax = cthetac
        endif
C
C Fetching yield stress from flow curve
        call ahard(sigmayield,hard,stateOld(i,1),props(33),nvalue)
C
C Eq. (11) - Calculating radius of elastic zone
        radius = sigmayield * ( one - ceta*(eta-eta0) ) * (cthetas
     1  + ((cthetaax-cthetas)* (gamma- ((gamma**(om+one))/(om+one)) ) ) )
     1  * facT * facStrrate * facDctlDmgVol
C
C Eq. (15) - Checking if yielding. The yielding factor facyld is set to be zero
C in the elastic zone, one if yielding
        facyld = zero
        if( q - radius .ge. zero ) facyld = one
C
C Eq. (20) - Calculating the tensor of the normal direction of plastic flow
C       Avoiding a divide by zero by avoiding that theta = zero
        if (theta.eq.zero) theta = theta + tol
C
        fac1 = con1 * three/(two*q)
        fac2 = facT * facStrrate * facDctlDmgVol
     1       * con1 * sigmayield*ceta*(cthetas+((cthetaax-cthetas)*(gamma
     1  -( (gamma**(om+one)) / (om+one) )))) * three*eta/(two*(q**2))
        fac3 = facT * facStrrate * facDctlDmgVol
     1       * con1 * sigmayield * (one-ceta*(eta-eta0)) * (cthetaax-cthetas)
     1  * (one -(gamma**om)) * (three*sqrt(three)/(two-sqrt(three)))
     1  * (tan(theta-pi/6.)/cos(theta-pi/6.)) * (one/(q*sin(three*theta)))
C
        on1 = fac1*ds1
     1  - fac2 * ds1
     1  - fac3 * (one/three+cos(three*theta)*ds1/(two*q)-three*((ds1
     1  **2)+(sig4**2)+(sig6**2))/(two*(q**2)))
C
        on2 = fac1*ds2
     1  - fac2 * ds2
     1  - fac3 * (one/three+cos(three*theta)*ds2/(two*q)-three*((sig4
     1  **2)+(ds2**2)+(sig5**2))/(two*(q**2)))
C
        on3 = fac1*ds3
     1  - fac2 * ds3
     1  - fac3 * (one/three+cos(three*theta)*ds3/(two*q)-three*((sig6
     1  **2)+(sig5**2)+(ds3**2))/(two*(q**2)))
C
        on4 = fac1*sig4
     1  - fac2 * sig4
     1  - fac3 * (cos(three*theta)*sig4/(two*q)-three*((ds1*sig4)
     1  +(ds2*sig4)+(sig6*sig5))/(two*(q**2)))
C
        on5 = fac1*sig5
     1  - fac2 * sig5
     1  - fac3 * (cos(three*theta)*sig5/(two*q)-three*((sig4*sig6)
     1  +(ds2*sig5)+(ds3*sig5))/(two*(q**2)))
C
        on6 = fac1*sig6
     1  - fac2 * sig6
     1  - fac3 * (cos(three*theta)*sig6/(two*q)-three*((ds1*sig6)
     1  +(sig4*sig5)+(ds3*sig6))/(two*(q**2)))
C
C Calculating trial yield function
        phitrial   = facyld * con1 * (q - radius) / facDctlDmgVol
C
C Calculating equivalent plastic strain
        deqps  = con1 * phitrial / (twomu + twoThirds * hard)
        stateNew(i,17) = deqps
C
C Updating equivalent plastic strain
        stateNew(i,1) = stateOld(i,1) + deqps
C
C Updating plastic strain (states 11-16)
        stateNew(i,11) = stateOld(i,11) + deqps * on1 / con1
        stateNew(i,12) = stateOld(i,12) + deqps * on2 / con1
        stateNew(i,13) = stateOld(i,13) + deqps * on3 / con1
        stateNew(i,14) = stateOld(i,14) + deqps * on4 / con1
        stateNew(i,15) = stateOld(i,15) + deqps * on5 / con1
        stateNew(i,16) = stateOld(i,16) + deqps * on6 / con1
C
C Updating equivalent plastic strain rate
        stateNew(i,22) = stateOld(i,22) + one
        stateNew(i,23) = stateOld(i,23) + deqps
        stateNew(i,24) = stateOld(i,24) + dt
        if (strrInc.gt.half .and. strrInc.lt.(one+half)) then
          stateNew(i,21) = deqps/dt
          stateNew(i,27) = deqps
        else
          stateNew(i,21) = (stateNew(i,23)/stateNew(i,24)
     1                   +  stateOld(i,21))/two
          stateNew(i,27) = (stateNew(i,23)/stateNew(i,22)
     1                   +  stateOld(i,27))/two
          if (stateNew(i,22).gt.strrInc-half .and. stateNew(i,22)
     1    .lt.strrInc+half) then
            stateNew(i,22) = zero
            stateNew(i,23) = zero
            stateNew(i,24) = zero
          endif
        endif
        stateNew(i,20) = deqps / dt
C
C Updating temperature
        if (tempOld(i).le.tol .and. tempOld(i).ge.-tol) then
          tempNew(i) = t0
     1    +(eta2/(density(i)*cp))*radius*deqps/facDctlDmgDev
        else
          tempNew(i) = tempOld(i)
     1    +(eta2/(density(i)*cp))*radius*deqps/facDctlDmgDev
        endif
C
C Calculating equivalent plastic strain to failure
        epsilonf = (half * ( d1*exp(zero-d2*eta) + d5*exp(zero-d6*eta))
     1                     - d3*exp(zero-d4*eta) )
     1                     * thetabar**2
     1           +  half * ( d1*exp(zero-d2*eta) - d5*exp(zero-d6*eta))
     1                     * thetabar
     1           +  d3*exp(zero-d4*eta)
C
C Ductile damage
        if (stateNew(i,1).ge.epsilonf .and. stateNew(i,1).gt.tol) then
          if (stateOld(i,3).le.tol .and. stateOld(i,3).ge.-tol) then
C           Registering the yield stress at the onset of damage for the first time
            stateNew(i,3) = radius
          else
C           Keeping the yield stress at the onset of damage for next step
            stateNew(i,3) = stateOld(i,3)
          endif
C         Updating damage fraction variable
          stateNew(i,2) = stateOld(i,2) + (stateNew(i,3)
     1                    /gf) * deqps
          stateNew(i,4) = two
        else
C
C         In case no damage is reached, yield stress at damage onset
C         and damage variable are kept for next step
          stateNew(i,3) = stateOld(i,3)
          stateNew(i,2) = stateOld(i,2)
        endif
C
C       Previous ductile damage
        if (stateOld(i,4).gt.(one+half)
     1  .and. stateOld(i,4).lt.(two+half)) then
          stateNew(i,4) = two
        endif
C
C Update stress
        fac1 = facDctlDmgVol * akappa * epsilontrace
        fac2 = facDctlDmgDev * twomu * deqps / con1
        sig1 = facDctlDmgDev * twomu * (stateNew(i,5)
     1       - epsilontrace/three - stateOld(i,11))
        sig2 = facDctlDmgDev * twomu * (stateNew(i,6)
     1       - epsilontrace/three -stateOld(i,12))
        sig3 = facDctlDmgDev * twomu * (stateNew(i,7)
     1       - epsilontrace/three -stateOld(i,13))
        sig4 = facDctlDmgDev * twomu * (stateNew(i,8) -stateOld(i,14))
        sig5 = facDctlDmgDev * twomu * (stateNew(i,9) -stateOld(i,15))
        sig6 = facDctlDmgDev * twomu * (stateNew(i,10)-stateOld(i,16))
        sig1 = sig1 + fac1 - fac2 * on1
        sig2 = sig2 + fac1 - fac2 * on2
        sig3 = sig3 + fac1 - fac2 * on3
        sig4 = sig4        - fac2 * on4
        sig5 = sig5        - fac2 * on5
        sig6 = sig6        - fac2 * on6
C       Calculating invariants of stress tensor
        SI1 =   sig1 + sig2 + sig3
        SI2 =   sig1*sig2-sig4*sig4
     1        + sig1*sig3-sig6*sig6
     1        + sig2*sig3-sig5*sig5
        SI3 =   sig1*(sig2*sig3-sig5*sig5)
     1        - sig4*(sig4*sig3-sig5*sig6)
     1        + sig6*(sig4*sig5-sig2*sig6)
C       Preparing subvalues for calculating the principal stresses values
        cosine2 = (two*SI1*SI1*SI1-three*three*SI1*SI2+three*three*three*SI3)
     1            /(two*(SI1*SI1-three*SI2)**(threehalves))
C       Assuring that -1<cosine2<1
        if (cosine2.gt.one) cosine2=one
        if (cosine2.lt.(zero-one)) cosine2=zero-one
        alpha2 = acos(cosine2)
C       Calculating the principal stress values
        SP1 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(alpha2/three)
        SP2 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(alpha2/three+twoThirds*pi)
        SP3 = SI1/three + twoThirds*sqrt(SI1*SI1-three*SI2)
     1                   *cos(twoThirds*pi-alpha2/three)
C       Fetching the highest of the principal stress values
        sigmamax = max(abs(SP1),abs(SP2),abs(SP3))
        stateNew(i,26) = sigmamax
C       Calculating new softening correction terms due to ductile damage
        stateNew(i,28)= one
        if (stateNew(i,2).ge.Dcrit) then
          facDctlDmgDev = zero
          facDctlDmgVol = zero
          stateNew(i,2) = Dcrit
          stateNew(i,28)= zero
        else
          facDctlDmgDev = one - stateNew(i,2)
          facDctlDmgVol = facDctlDmgDev
        endif
C       Elements under hydrostatic compression don't experience spherical damage
        if (sig1+sig2+sig3.lt.zero .and. stateNew(i,2).lt.Dcrit) then
          facDctlDmgVol = one
        endif
        fac1 = facDctlDmgVol * akappa * epsilontrace
        fac2 = facDctlDmgDev * twomu * deqps / con1
        sig1 = facDctlDmgDev * twomu * (stateNew(i,5)
     1       - epsilontrace/three - stateOld(i,11))
        sig2 = facDctlDmgDev * twomu * (stateNew(i,6)
     1       - epsilontrace/three -stateOld(i,12))
        sig3 = facDctlDmgDev * twomu * (stateNew(i,7)
     1       - epsilontrace/three -stateOld(i,13))
        sig4 = facDctlDmgDev * twomu * (stateNew(i,8) -stateOld(i,14))
        sig5 = facDctlDmgDev * twomu * (stateNew(i,9) -stateOld(i,15))
        sig6 = facDctlDmgDev * twomu * (stateNew(i,10)-stateOld(i,16))
        stressNew(i,1) = sig1 + fac1 - fac2 * on1
        stressNew(i,2) = sig2 + fac1 - fac2 * on2
        stressNew(i,3) = sig3 + fac1 - fac2 * on3
        stressNew(i,4) = sig4        - fac2 * on4
        stressNew(i,5) = sig5        - fac2 * on5
        stressNew(i,6) = sig6        - fac2 * on6
C       Brittle damage
        if (sigmamax.ge.sigdmg .and. stressNew(i,1).gt.zero) then
          stressNew(i,1) = zero
          stressNew(i,2) = zero
          stressNew(i,3) = zero
          stressNew(i,4) = zero
          stressNew(i,5) = zero
          stressNew(i,6) = zero
          stateNew(i,4)  = one
          stateNew(i,28) = zero
        endif
C
C
  100 continue
C
      return
      end
C
C
      subroutine ahard(sigmayield,hard,eqplas,table,nvalue)
C
      include 'vaba_param.inc'
      dimension table(2,nvalue)
C
C     Set yield stress to last value of table, hardening to zero
      sigmayield=table(1,nvalue)
      hard=0.0
C
C     If more than one entry, search table
      if(nvalue.gt.1) then
        do 10 k1=1,nvalue-1
          eqpl1=table(2,k1+1)
          if(eqplas.lt.eqpl1) then
            eqpl0=table(2,k1)
            if(eqpl1.le.eqpl0) then
              write(6,7)
 7            format(//,30X,'***ERROR - PLASTIC STRAIN MUST BE ',
     1               'ENTERED IN ASCENDING ORDER,')
C
C             Subroutine XIT terminates execution and closes all files
C              call XIT
            endif
            deqpl=eqpl1-eqpl0
            sigmayield0=table(1,k1)
            sigmayield1=table(1,k1+1)
            dsigmayield=sigmayield1-sigmayield0
            hard=dsigmayield/deqpl
            sigmayield=sigmayield0+(eqplas-eqpl0)*hard
            goto 20
          endif
 10     continue
 20     continue
        if(eqplas.gt.table(2,nvalue)) then
          hard=(table(1,nvalue)-table(1,nvalue-1))
     1        /(table(2,nvalue)-table(2,nvalue-1))
          sigmayield=table(1,nvalue)+(eqplas-table(2,nvalue))*hard
        endif
      endif
      return
C
C Iteration ends here
      end
