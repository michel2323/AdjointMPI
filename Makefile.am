SUBDIRS = src 

EXTRA_DIST = examples interfaces

dist-hook:
		rm -rf `find $(distdir)/examples -type d -name .svn`
		rm -rf `find $(distdir)/interfaces -type d -name .svn`

