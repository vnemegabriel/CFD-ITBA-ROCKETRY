DOC = aconcagua-user-guide
TEX = lualatex -interaction=nonstopmode -halt-on-error

$(DOC).pdf: $(DOC).tex settings.tex references.bib $(wildcard chapters/*.tex) $(wildcard figures/pdf/*.pdf)
	$(TEX) $(DOC) > /dev/null || { tail -40 $(DOC).log; exit 1; }
	bibtex $(DOC) > /dev/null || { cat $(DOC).blg; exit 1; }
	$(TEX) $(DOC) > /dev/null
	$(TEX) $(DOC) > /dev/null || { tail -40 $(DOC).log; exit 1; }
	@grep -E 'Warning.*(undefined|multiply)' $(DOC).log || true

figures:
	node figures/render.mjs

clean:
	rm -f $(DOC).{aux,bbl,blg,log,out,toc,lof,lot}

.PHONY: figures clean
