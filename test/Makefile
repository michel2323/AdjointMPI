TESTS=$(shell find . -type d)

TEST_DIR=.

check: 
	for i in $(TESTS); do if [ $$i != "." ]; then cd ${TEST_DIR}/$$i && echo -n $$i && make check && cd ..; fi; done

clean: 
	for i in $(TESTS); do if [ $$i != "." ]; then cd ${TEST_DIR}/$$i && make clean && cd ..; fi; done
