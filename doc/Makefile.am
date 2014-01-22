if HAVE_DOXYGEN

doxyfile.stamp:
	@echo "Generating documentation using doxygen..."
	doxygen Doxyfile.in >& doxygen.log
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp

all: doxyfile.stamp
clean:
	rm -rf ./latex
	rm -rf ./html
	rm -f doxyfile.stamp
endif 
