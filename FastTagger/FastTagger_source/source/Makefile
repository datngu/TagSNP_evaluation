CC = g++
CFLAGS = -std=c++11 -Wall
OBJECTS = main.o gen_rules.o loaddata.o rule_io.o select_tag.o verify.o
TARGET = FastTaggerV2

$(TARGET) : $(OBJECTS) Makefile
	$(CC) $(CFLAGS) $(OBJECTS) -o $(TARGET)

gen_rules.o: gen_rules.cpp global.h
	$(CC) -c $(CFLAGS) gen_rules.cpp

loaddata.o: loaddata.cpp global.h
	$(CC) -c $(CFLAGS) loaddata.cpp

main.o: main.cpp global.h
	$(CC) -c $(CFLAGS) main.cpp

rule_io.o: rule_io.cpp global.h
	$(CC) -c $(CFLAGS) rule_io.cpp

select_tag.o: select_tag.cpp global.h
	$(CC) -c $(CFLAGS) select_tag.cpp

verify.o: verify.cpp global.h
	$(CC) -c $(CFLAGS) verify.cpp

clean :
	rm *.o
