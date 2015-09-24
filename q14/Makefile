HOSTS = antares hyades praesepe virgo astronomy astrodisk

default:  link html tex clean

link:
	ln -fs $(CURDIR)/resources ~/bin
	ln -fs $(CURDIR)/resources_html ~/bin
	ln -fs $(CURDIR)/resources_tex ~/bin

html:
	rm resources.html 
	getresources_html $(HOSTS) >resources.html

tex:
	rm resources.tex
	getresources_tex $(HOSTS) >resources.tex
	pdflatex resources
	-rm resources.aux resources.log

clean:
	-rm ~/bin/resources ~/bin/resources_html ~/bin/resources_tex
