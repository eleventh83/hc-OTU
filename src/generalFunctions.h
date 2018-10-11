#ifndef GENERALFUNCTIONS_H_
#define GENERALFUNCTIONS_H_


std::string& itostr(int n, std::string& s);

std::string int2str(int n);

std::string dbl2str(double n);

std::string replace(const std::string old_string, const std::string new_string, std::string target_string);

std::string &ltrim(std::string &s);
std::string &rtrim(std::string &s);
std::string &bothtrim(std::string &s);


std::string trim(std::string str);

std::string toUpperCase(std::string str);

std::vector<std::string> StringSplit(std::string strOrigin, std::string strTok, int limit);


#endif
