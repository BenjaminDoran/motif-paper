usage="$(basename "$0") [-h] [outputFile] -- Builds Motif Paper"
while getopts ':hs:' option; do
    case "$option" in
        h) echo "$usage"
           exit
           ;;
    esac
done


# input files
read -r -d '' INPUT << EOM
./sections/paper-titlepage.md 
./sections/paper-abstract.md 
./sections/paper-intro.md 
./sections/paper-background.md 
./sections/paper-methodology.md 
./sections/paper-evaluation.md 
./sections/paper-conclusion.md 
./sections/paper-references.md
EOM

# output file
[[ $1 ]] && OUTPUT=$1 || OUTPUT="BenjaminDoran_MotifPaper.docx"

# flags and style
read -r -d '' FLAGS << EOM
--smart --standalone 
-V geometry:margin=1in 
--reference-docx=./template/custom-reference.docx
EOM

# Citation settings
read -r -d '' REFERENCE << EOM
--bibliography ./reference/Motif.bib 
--filter pandoc-citeproc 
--csl=./reference/ieee.csl
EOM

# make paper
pandoc $REFERENCE $INPUT $FLAGS -o $OUTPUT