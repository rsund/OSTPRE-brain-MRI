#! /bin/bash
# This takes subjses file as an input and generates the required directories and files to run synthseg segmentation 
if [ "$1" = "" ]; then
  echo usage:  $0 bids_directory_sub_ses
  exit
fi
sessionid=$( basename $1 ) # get session
pathnametmp=$(dirname $1 )
subjid=$( basename ${pathnametmp}) # get subjectid
pathname=$(dirname ${pathnametmp})
subjsession="${subjid}_${sessionid}" 
t1=`ls $1/anat/sub*_T1w.nii.gz 2>/dev/null`  # this is complete filename 
t1file=$( basename ${t1} )
# remember to add a 3T option if that seems like a valid consideration 
bids_folder_cross="${pathname}/derivatives/synthseg2.0/${subjid}"     # define BIDS path for cross-sectional data
bids_folder_ses="${bids_folder_cross}/${sessionid}"
if [ ! -d ${bids_folder_cross} ]; then   # make the directory if it does not exist
    mkdir ${bids_folder_cross}
fi     
if [ ! -d ${bids_folder_ses} ]; then   # make the directory if it does not exist
    mkdir ${bids_folder_ses}
fi 
t1filename="${pathname}/synthseg_t1file.txt"
segfilename="${pathname}/synthseg_segfile.txt"
csvfilename="${pathname}/synthseg_csvfile.txt"
qcfilename="${pathname}/synthseg_qcfile.txt"
resamplefilename="${pathname}/synthseg_resamplefile.txt"

echo "${t1}" >> ${t1filename}
echo "${bids_folder_ses}/synthseg_${t1file}" >> ${segfilename}
echo "${bids_folder_ses}/synthseg-vol_${subjsession}.csv" >> ${csvfilename} 
echo "${bids_folder_ses}/synthseg-qc_${subjsession}.csv" >> ${qcfilename}
echo "${bids_folder_ses}/resample_${t1file}" >> ${resamplefilename}


# mri_synthseg --i ${t1} --o ${segfile} --vol ${csvfile}  --qc ${qcfile} --resample ${resamplefile} --cpu

