default: test

montgomery.a: montgomery.o
	ar rcs $@ $<

test: test.o montgomery.a

clean:
	$(RM) -f *.o *.a test
