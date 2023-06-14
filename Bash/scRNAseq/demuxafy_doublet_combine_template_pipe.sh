#!/bin/bash

# Running Demuxafy
echo $'\nRunning Demuxafy: Doublet Detecting & Combine_Results . . . \n'

#  Input Paths (CHANGE HERE FOR EACH SAMPLE)
#  ===================================================================
COUNTS=

# Output Directories
FREEMUXLET_OUTDIR=
SOUPORCELL_OUTDIR=
VIREO_OUTDIR=
SCSPLIT_OUTDIR=

SCDBLFINDER_OUTDIR=
SCDS_OUTDIR=
Combine_Results_OUTDIR=
#  ===================================================================

# Print input paths
echo "----------------------------------------"
echo $'\n'
echo $' COUNTS: $COUNTS\n'
echo $' FREEMUXLET_OUTDIR: $FREEMUXLET_OUTDIR\n'
echo $' SOUPORCELL_OUTDIR: $SOUPORCELL_OUTDIR\n'
echo $' VIREO_OUTDIR: $VIREO_OUTDIR\n'
echo $' SCSPLIT_OUTDIR: $SCSPLIT_OUTDIR\n'
echo $' SCDBLFINDER_OUTDIR: $SCDBLFINDER_OUTDIR\n'
echo $' SCDS_OUTDIR: $SCDS_OUTDIR\n'
echo $' Combine_Results_OUTDIR: $Combine_Results_OUTDIR\n'
echo $'\n'
echo "----------------------------------------"


# calculte run time
start=$(date +%s)

#  ===================================================================
#  SCDBLFINDER
#  ===================================================================
echo $'\nRunning SCDBLFINDER . . . \n'

# calculte run time
start_SCDBLFINDER=$(date +%s)

singularity exec Demuxafy.sif SCDBLFINDER.R \
    -o $SCDBLFINDER_OUTDIR \
    -t $COUNTS


echo $'\SCDBLFINDER complete! \n'

end_SCDBLFINDER=$(date +%s)
runtime_SCDBLFINDER=$((end_SCDBLFINDER-start_SCDBLFINDER))

#  ===================================================================
#  SCDS
#  ===================================================================
echo $'\nRunning SCDS . . . \n'

# calculte run time
start_SCDS=$(date +%s)

singularity exec Demuxafy.sif scds.R \
    -o $SCDS_OUTDIR \
    -t $COUNTS


echo $'\nSCDS complete! \n'

end_SCDS=$(date +%s)
runtime_SCDS=$((end_SCDS-start_SCDS))

#  ===================================================================
#  Combine_Results
#  ===================================================================
echo $'\nRunning Combine_Results . . . \n'

# calculte run time
start_Combine_Results=$(date +%s)


singularity exec Demuxafy.sif Combine_Results.R \
  -o $Combine_Results_OUTDIR/combined_results.tsv \
  --FREEMUXLET $FREEMUXLET_OUTDIR \
  --souporcell $SOUPORCELL_OUTDIR \
  --VIREO $VIREO_OUTDIR \
  --SCSPLIT $SCSPLIT_OUTDIR \
  --SCDBLFINDER $SCDBLFINDER_OUTDIR \
  --SCDS $SCDS_OUTDIR \
  --method "MajoritySinglet"

echo $'\nCombine_Results complete! \n'

end_Combine_Results=$(date +%s)
runtime_Combine_Results=$((end_Combine_Results-start_Combine_Results))

# total runtime
end = $(date +%s)
runtime = $((end-start))

# ===================================================================

# Print runtimes
echo "----------------------------------------"
echo $'Pipeline runtimes: \n'
echo "----------------------------------------"
echo $'Demuxafy: $runtime_demuxafy seconds \n'
echo $'SCDBLFINDER: $runtime_SCDBLFINDER seconds \n'
echo $'SCDS: $runtime_SCDS seconds \n'
echo $'Combine_Results: $runtime_Combine_Results seconds \n'
echo "----------------------------------------"
echo $'Total runtime: $runtime seconds \n'
echo "----------------------------------------"
