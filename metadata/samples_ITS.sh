#!/bin/bash
# Runs in bash on Mac OS X 13.0.1

# Reformat ITS metadata from ANL and fix error

# File names
INPUT="220520_Doherty_ITS_Closeout_220518_MODIFIED.txt"
TEMP1="samples_ITS.temp1"
TEMP2="samples_ITS.temp2"
TEMP3="samples_ITS.temp3"
OUTPUT="samples_ITS.tsv"

# Input file appears to have iso-8859-1 character encoding
#     file -I ${INPUT}

# Discard non UTF-8 characters
iconv -f utf-8 -t utf-8 -c ${INPUT} >${TEMP1}

# Alternatively, replace Mac encoded characters with UTF-8.
# This retains the degree symbol but might well prove problematic later.
#     iconv -f MAC -t utf-8 -c ${INPUT} >${TEMP1}

# Clean up line breaks and extra tabs
tr '\r' '\n' <${TEMP1} | sed 's/\t\t//' >${TEMP2}

# Select metadata to keep
# Removes samples from separate respiration study (unrelated to SERDP RDX project
# but sequenced together).
# grep -v ".*C-.*kpa-R.*\t" ${TEMP2} ${TEMP3}

# Keep respiration samples in for now to examine how many reads go to those samples
cp ${TEMP2} ${TEMP3}

# Process remaining file contents
awk 'BEGIN {FS = "\t"; OFS="\t"}
NR==1 {
    print "SampleID", "Dataset", "Project", "BarcodeSequence", "LinkerPrimerSequence", "BarcodePlate", "Well"
}
!(NR==1) && !($1=="") {
    sample = $1;
    dataset = "ITS";
    if (sample ~ /.*C-.*kpa-R.*/)
    {
        gsub(/-/, "_", sample);
        project = "respiration";
    }
    else if (sample ~ /^SERDP EB/)
    {
        sub(/ /, "_", sample);
        project = "serdp_blanks";
    }
    else
    {
        gsub(/\+/, "y", sample); gsub(/-/, "n", sample);
        project = "serdp";
    }
    print sample, dataset, project, $2, $3, $4, $5
}' ${TEMP3} >${OUTPUT}

# Clean up
rm ${TEMP1} ${TEMP2} ${TEMP3}
unset INPUT TEMP1 TEMP2 TEMP3 OUTPUT
