#CC = /opt/intel/bin/icc
CC = g++
#OMPFLAGS = -openmp
OMPFLAGS = -fopenmp
CCFLAGS = -O3 -g -fPIC -Wall
CONSTFLAGS = -Wno-write-strings

TARGET = hc-OTU
OBJECTS = ssw.o ssw_cpp.o AlignNoGap.o AlignNoGapWrap.o AssembleGlobal.o CLHeader.o CLSQFastaReader.o ContigCluster.o Contig.o FastaSeq.o Fasta.o Fastq.o CLSQSequence.o generalFunctions.o KmerCluster.o KmerSimple.o KmerDist.o KmerHitLine.o needle.o Nucleotide.o Seqtools.o TBCclass.o TKmerCluster.o Tools.o util.o



all: main.cpp $(OBJECTS)
	$(CC) $(CCFLAGS) $(OMPFLAGS) $(CONSTFLAGS) $< $(OBJECTS) -o $(TARGET)

ssw.o: ssw.c
	$(CC) $(CCFLAGS) -c $<

ssw_cpp.o: ssw_cpp.cpp
	$(CC) $(CCFLAGS) -c $<

AlignNoGap.o: AlignNoGap.cpp
	$(CC) $(CCFLAGS) -c $<

AlignNoGapWrap.o: AlignNoGapWrap.cpp
	$(CC) $(CCFLAGS) -c $<

AssembleGlobal.o: AssembleGlobal.cpp
	$(CC) $(CCFLAGS) -c $<

CLHeader.o: CLHeader.cpp
	$(CC) $(CCFLAGS) -c $<

CLSQFastaReader.o: CLSQFastaReader.cpp
	$(CC) $(CCFLAGS) -c $<

ContigCluster.o: ContigCluster.cpp
	$(CC) $(CCFLAGS) -c $<

Contig.o: Contig.cpp
	$(CC) $(CCFLAGS) -c $<

Fasta.o: Fasta.cpp
	$(CC) $(CCFLAGS) -c $<

Fastq.o: Fastq.cpp
	$(CC) $(CCFLAGS) -c $< 

CLSQSequence.o: CLSQSequence.cpp
	$(CC) $(CCFLAGS) -c $< 

FastaSeq.o: FastaSeq.cpp
	$(CC) $(CCFLAGS) -c $<

generalFunctions.o: generalFunctions.cpp
	$(CC) $(CCFLAGS) -c $<

KmerCluster.o: KmerCluster.cpp
	$(CC) $(OMPFLAGS) $(CCFLAGS) -c $<

KmerSimple.o: KmerSimple.cpp
	$(CC) $(OMPFLAGS) $(CCFLAGS) -c $<

KmerDist.o: KmerDist.cpp
	$(CC) $(CCFLAGS) -c $<

KmerHitLine.o: KmerHitLine.cpp
	$(CC) $(CCFLAGS) -c $<

needle.o: needle.cpp
	$(CC) $(CCFLAGS) -c $<

Nucleotide.o: Nucleotide.cpp
	$(CC) $(CCFLAGS) -c $<

Seqtools.o: Seqtools.cpp
	$(CC) $(CCFLAGS) -c $<

TBCclass.o: TBCclass.cpp
	$(CC) $(OMPFLAGS) $(CCFLAGS) -c $<

TKmerCluster.o: TKmerCluster.cpp
	$(CC) $(CCFLAGS) -c $<

Tools.o: Tools.cpp
	$(CC) $(CCFLAGS) -c $<

util.o: util.cpp
	$(CC) $(CCFLAGS) -c $<


clean:
	rm -f hc-OTU *.o


