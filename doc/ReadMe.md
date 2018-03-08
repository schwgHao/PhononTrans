## phonon transport
This is the phonon transport part in SMEAGOL. 
--- Mar. 8, 2018 ---
--- by Hao Wang  ---

### compile
The DFT part of SMEAGOL was written with Fortran. To reduce the coupling to the initial code and better encapsule , C++ was used to encode this part. Some characteristic of Fortran2003 and C++11 were used, especially in building the interface. The following compile labels are necessary in arch.make,
- CXX = g++ #(newer than 4.8)
- CCFLAGS = -Wall -Wextra -std=c++11
- LDFLAGS = $(FFLAGS) libstdc++.so.6.0.17 #(must newer than 17)

### execute
- run script "run_phononTrans.sh slabel.fdf" to construct a huge non-peiordic super cell for lead
    - no further label is needed, the script will add them automatically
    - systemlabel.LeadInfo will be generate to contain calculating information
    - slabelFC.fdf then can be used as input file for smeagol
- run "pathforsmeagol/smeagol < slabelFC.fdf >slabelFc.out" for lead calculations
    - bulksgf.lft, bulksgf.rt will be generated. Copy them to the transport calculation directory.

- in transport directory, 
    - MD.FCfirst   14
    - MD.FClast    28
    should be added to boundary the central scattering region,
- run "run_phononTrans.sh slabel.fdf E" to construct supercell for transport calculation
    - generate systemlabel.EMInfo, slabelFC.fdf
- run "pathforsmeagol/smeagol < slabelFC.fdf >slabelFc.out" for transport calculations
    - generate:
        - systemlabel.PhTrc  // transmission
        - systemlabel.THERMPH  // thermal conductance
        - systemlabel.BulkDOS  // system DOS before and after connecting molecules to leads
        - systemlabel.EigenMode // eigenmode and eigenvectors of local vibration
        - systemlabel.EigenChannel // Eigen conduction channel
- calulate decay rate of local vibration
    - copy systemlabel.EigenMode to lead directory and add 
        vibDecayRate        T
      to slabelFc.fdf and run smeagol
    - copy bulksgf.lft, bulksgf.rt to transport directory and add
        vibDecayRate        T
      to slabelFc.fdf in transport directory and run smeagol
    - generate
        - systemlabel.EigenModeRate  // decay rate of local vibration
- also, the demanding FC matrix calculations could also be performed with siesta, the label
    - MD.FCRead           T
    - MD.FCAtomRestart    10000
    should be added into the fdf file.

### end


