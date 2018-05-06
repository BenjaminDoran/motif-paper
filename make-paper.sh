usage="$(basename "$0") [-h] [-o outputFile] -- Builds Motif Paper"
while getopts ':hdos:' option; do
    case "$option" in
        h)  echo "$usage"
            exit
        ;;
        d) debug=true ;;
        o)  shift
            OUTPUT=$1
        ;;
    esac
done


# input files
read -r -d '' INPUT << EOM
./sections/paper-titlepage.md 
./sections/paper-abstract.md 
./sections/paper-intro.md 
./sections/paper-dip.md 
./sections/paper-methodology.md 
./sections/paper-evaluation.md 
./sections/paper-conclusion.md 
./sections/paper-references.md
EOM

# output file
[[ $OUTPUT ]] || OUTPUT="BenjaminDoranMotifPaper.docx"

# flags and style
read -r -d '' FLAGS << EOM
--standalone -f markdown -V geometry:margin=1in 
--mathjax
--css=./template/custom.css
--reference-doc=./template/custom-reference
EOM
[[ $(echo $OUTPUT | grep -o '[^\.]\+$') == 'tex' ]] && FLAGS=$FLAGS.tex
[[ $(echo $OUTPUT | grep -o '[^\.]\+$') == 'docx' ]] && FLAGS=$FLAGS.docx


# Citation settings
read -r -d '' REFERENCE << EOM
--bibliography ./reference/Motif.bib 
--filter pandoc-citeproc 
--csl=./reference/ieee.csl
EOM

# make paper
[[ $debug ]] && echo "/usr/bin/pandoc $REFERENCE $INPUT $FLAGS -o $OUTPUT"
/usr/bin/pandoc $REFERENCE $INPUT $FLAGS -o $OUTPUT