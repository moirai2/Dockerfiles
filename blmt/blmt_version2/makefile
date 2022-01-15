all : clean faa2srt srt2lcp ngrams proteinCount proteinNGram 
CFLAGS  =  -O2
LDFLAGS =
CC = gcc

clean-all: 
	rm -rf *.o faa2srt srt2lcp ngrams proteinCount proteinNGram
clean:
	rm -rf *.o

mylib.o: 
	$(CC) $(CFLAGS) $(LDFLAGS) -c mylib.c
faa2srt.o: 
	$(CC) $(CFLAGS) $(LDFLAGS) -c faa2srt.c
srt2lcp.o:
	$(CC) $(CFLAGS) $(LDFLAGS) -c srt2lcp.c
ngrams.o: 
	$(CC) $(CFLAGS) $(LDFLAGS) -c ngrams.c
proteinCount.o:
	$(CC) $(CFLAGS) $(LDFLAGS) -c proteinCount.c
proteinNGram.o:
	$(CC) $(CFLAGS) $(LDFLAGS) -c proteinNGram.c

faa2srt: faa2srt.o mylib.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o faa2srt faa2srt.o mylib.o
srt2lcp: srt2lcp.o mylib.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o srt2lcp srt2lcp.o mylib.o
ngrams: ngrams.o mylib.o  	
	$(CC) $(CFLAGS) $(LDFLAGS) -o ngrams ngrams.o mylib.o
proteinCount: proteinCount.o mylib.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o proteinCount proteinCount.o mylib.o
proteinNGram: proteinNGram.o mylib.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o proteinNGram proteinNGram.o mylib.o
