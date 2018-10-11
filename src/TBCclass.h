#ifndef TBCCLASS_H_
#define TBCCLASS_H_


class TBCclass {
public:
	TBCclass();
	~TBCclass();
	void assembleWithRead(std::string, double , int, double);

	int nSingleton;
	int nOtus;
	std::string inmc;
	int nIter;
};

#endif /* TBCCLASS_H_ */
