TESTS=$(shell find . -type d)

TEST_DIR=.

check: 
	for i in `ls -l | egrep '^d'  | awk ' { print $$9 } '` ; do cd $$i && make check && cd .. ; done

clean: 
	for i in `ls -l | egrep '^d'  | awk ' { print $$9 } '` ; do cd $$i && make clean && cd .. ; done
