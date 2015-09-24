HOSTS = antares hyades praesepe virgo astronomy astrodisk

default:  link html tex clean

link:
	ln -fs $(CURDIR)/resources ~/bin

html:
	rm resources.html 
	getresources html $(HOSTS) >resources.html

tex:
	rm resources.tex
	getresources tex $(HOSTS) >resources.tex
	pdflatex resources
	-rm resources.aux resources.log

clean:
	-rm ~/bin/resources
