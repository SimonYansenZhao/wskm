APP=wskm
VER=1.4.0
DATE=$(shell date +%Y-%m-%d)

######################################################################
# R Package Management

.PHONY: check
check: clean build
	R CMD check --as-cran --check-subdirs=yes $(APP)_$(VER).tar.gz

.PHONY: build
build: 
	perl -pi -e 's|^Version: .*|Version: $(VER)|' package/DESCRIPTION
	perl -pi -e 's|^Date: .*|Date: $(DATE)|' package/DESCRIPTION
	R CMD build package

.PHONY: install
install: build
	R CMD INSTALL $(APP)_$(VER).tar.gz

.PHONY: rebuild
rebuild: build install

########################################################################
# Misc

.PHONY: clean
clean:
	@rm -vf ./*~
	@rm -vf package/*/*~
	@rm -vf package/*/.Rhistory
	@rm -rf package.Rcheck rattle.Rcheck
