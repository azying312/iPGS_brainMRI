#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
REPODIR=$(dirname $(dirname ${SRCDIR}))

source ${REPODIR}/paths.sh
analysis_name=$(basename ${SRCDIR})

ml load pigz

# output
#SQL_dump_f="${ukb21942_d}/pheno/tmp/${analysis_name}/tab_dump.tsv.gz"
SQL_dump_f="${ukb21942_user}/pheno/tmp/${analysis_name}/tab_dump.tsv.gz"
echo "SQL_dump_f: ${SQL_dump_f}"

# input
long_f="${ukb21942_d}/pheno/basket4043067/ukb672713.long.tsv.gz"
echo "long_f: ${long_f}"

fields_f="${SRCDIR}/misc/mri_traits_449_names.tsv"
#fields_f="${SRCDIR}/misc/short_fields.tsv"
#fields_f="${SRCDIR}/misc/one_fields.tsv"
echo "fields: ${fields_f}"

# main
if [ ! -d $(dirname "${SQL_dump_f}") ] ; then mkdir -p $(dirname "${SQL_dump_f}") ; fi

regex_pattern=$(cat ${fields_f} | awk '(NR>1){print "^" $1 "$" }' | tr '\n' '|' | rev | cut -c2- | rev )
echo "regex_pattern: ${regex_pattern}"

zcat ${long_f} \
| awk -v ptn=${regex_pattern} '((NR == 1) || ($2 ~ ptn))' \
| pigz -9 -p6 > ${SQL_dump_f}



# Check
if [ -s ${SQL_dump_f} ]; then
    echo "Output file ${SQL_dump_f} created successfully."
    echo "Number of lines in the output file: $(zcat ${SQL_dump_f} | wc -l)"
    echo "Size of the output file: $(du -h ${SQL_dump_f} | cut -f1)"
else
    echo "Error: Output file ${SQL_dump_f} is empty or was not created."
fi
