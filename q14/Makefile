HOSTS = antares hyades praesepe virgo astronomy astrodisk

default:  link html tex

link:
	ln -fs $(CURDIR)/resources ~/bin
	ln -fs $(CURDIR)/getresources_html ~/bin
	ln -fs $(CURDIR)/getresources_tex ~/bin

html:
	rm resources.html 
	getresources_html $(HOSTS) >resources.html

tex:
	rm resources.tex
	getresources_tex $(HOSTS) >resources.tex
	pdflatex resources
	-rm resources.aux resources.log

