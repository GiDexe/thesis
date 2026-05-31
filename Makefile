.PHONY: results paper rerun clean protect unprotect

paper: results
	latexmk -pdf -interaction=nonstopmode paper.tex

results:
	quarto render MAIN.qmd

rerun:
	RERUN=true Rscript estimation.r
	USE_RERUN=true quarto render MAIN.qmd
	latexmk -pdf -interaction=nonstopmode paper.tex

# Delete generated artifacts only — never output/ or output_lagged/
clean:
	rm -rf output_rerun
	rm -f  assets/tables/*.tex assets/figures/*.pdf assets/figures/*.png
	rm -f  results.pdf MAIN.pdf
	latexmk -C paper.tex 2>/dev/null || rm -f paper.pdf paper.aux paper.log paper.out
	rm -rf MAIN_files MAIN_cache .quarto

# Lock / unlock the reference estimates
protect:
	chmod -R a-w output output_lagged

unprotect:
	chmod -R u+w output output_lagged