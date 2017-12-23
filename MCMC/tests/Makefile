CC = g++

CPPFLAGS = -std=c++11

%.d: %.c
	@set -e \
	$(CC) -M $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\l.o $@ : ,g' < $@.$$$$ > $@; 
