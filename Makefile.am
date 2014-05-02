SUBDIRS = src 

EXTRA_DIST = examples interfaces test doc



doxyfile.stamp:
	@echo "Generating documentation using doxygen..."
	cd doc ; doxygen Doxyfile.in >& doxygen.log ; echo Timestamp > ../doxyfile.stamp

CLEANFILES = doxyfile.stamp

doc: doxyfile.stamp
clean-local:
	rm -rf ./doc/latex
	rm -rf ./doc/html
	rm -f doxyfile.stamp
