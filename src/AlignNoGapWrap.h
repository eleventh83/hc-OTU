#ifndef ALIGNNOGAPWRAP_H_
#define ALIGNNOGAPWRAP_H_

class AlignNoGapWrap:AlignNoGap{
public:
	AlignNoGapWrap();
	AlignNoGapWrap(std::string qseq, std::vector<CLSQSequence> seqsList, int startidx, int endidx, int th_overlap);
	~AlignNoGapWrap();

	std::vector<int> match_idx;
	std::vector<CLSQSequence> seqList;
	int idxstart;
	int idxend;

	void setQseq(std::string qseq);
	void setRange(int startidx, int endidx);
	void setMinOverlap(int th_overlap);
	void setSeqList(std::vector<CLSQSequence> seqsList);

	void clear();
	std::vector<int> call();
};

#endif /* ALIGNNOGAPWRAP_H_ */

