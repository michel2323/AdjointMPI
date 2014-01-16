if HAVE_DOXYGEN

doxyfile.stamp:
	doxygen Doxyfile.in
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp

all: doxyfile.stamp
clean:
	rm -rf ./latex
	rm -rf ./html
	rm -f doxyfile.stamp
endif 
