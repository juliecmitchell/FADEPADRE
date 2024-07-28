# toplevel Makefile 

###########################################

DIST	= fade

###########################################

MAKE = make

default: install 

build:
	@echo Making executable in src...
	@cd src; ../conf/configure -checkconfig;
	cd src; $(MAKE)

install:
	@echo Installing executable in src...
	@cd src; if [ ! -f Makefile ]; then ../conf/configure -checkconfig; fi;
	cd src; $(MAKE) install

clean:
	@echo Cleaning in src...
	@cd src; ../conf/configure -checkconfig; 
	cd src; $(MAKE) clean

realclean:
	@echo Making realclean in src...
	@cd src; ../conf/configure -checkconfig;
	cd src; $(MAKE) realclean

# Build distribution archive
dist:
	tar cvf $(DIST).tar ./Makefile ./README ./conf ./doc ./src ./bin
	/bin/rm -f $(DIST).tar.Z
	compress $(DIST).tar

