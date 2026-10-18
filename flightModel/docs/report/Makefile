DOC = aconcagua-cfd-report
TEX = lualatex -interaction=nonstopmode -halt-on-error

$(DOC).pdf: $(DOC).tex ../guide/references.bib ../../flightData/cdFlight.png
	$(TEX) $(DOC) > /dev/null || { tail -40 $(DOC).log; exit 1; }
	bibtex $(DOC) > /dev/null || { cat $(DOC).blg; exit 1; }
	$(TEX) $(DOC) > /dev/null
	$(TEX) $(DOC) > /dev/null || { tail -40 $(DOC).log; exit 1; }
	@grep -E 'Warning.*(undefined|multiply)' $(DOC).log || true

clean:
	rm -f $(DOC).{aux,bbl,blg,log,out}

.PHONY: clean
