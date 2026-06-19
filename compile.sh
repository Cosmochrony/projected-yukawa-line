#!/bin/bash
# Projected Yukawa Line note compilation script
set -e
TEX_FILE="tex/ProjectedYukawaLine.tex"
OUTPUT_DIR="out"
MAIN_NAME="ProjectedYukawaLine"
mkdir -p "$OUTPUT_DIR"
export TEXINPUTS=".:./tex:${TEXINPUTS}"
pdflatex -file-line-error -interaction=nonstopmode -output-directory="$OUTPUT_DIR" "$TEX_FILE"
( cd "$OUTPUT_DIR" && BSTINPUTS="../tex:${BSTINPUTS}" BIBINPUTS="../tex:${BIBINPUTS}" bibtex "$MAIN_NAME" ) || true
pdflatex -file-line-error -interaction=nonstopmode -output-directory="$OUTPUT_DIR" "$TEX_FILE"
pdflatex -file-line-error -interaction=nonstopmode -output-directory="$OUTPUT_DIR" "$TEX_FILE"
echo "Done: $OUTPUT_DIR/$MAIN_NAME.pdf"
