#ifndef TOOLS_H_
#define TOOLS_H_

namespace Tools{

	const std::string alphabet = "ABCDEFGHIJGLMNOPQRSTUVWXYZ1234567890";


	double calculate(std::string s1, std::string s2);
	std::string padEndGaps(std::string seq, char fill);
	std::vector<Contig> ContigFile2ContigList(std::string contig_fn);
	std::string padString(std::string s, int n, char c, bool paddingLeft);
	double getTimeStampNow();
}

#endif /* TOOLS_H_ */
