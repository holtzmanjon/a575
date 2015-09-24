HOSTS = antares hyades praesepe virgo astronomy astrodisk

default:  html tex

html:
	rm resources.html 
	getresources $(HOSTS) >resources.html

tex:
	rm resources.tex
	getresources_tex $(HOSTS) >resources.tex
	pdflatex resources

