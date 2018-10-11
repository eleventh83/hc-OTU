#include <algorithm>
#include <sstream>
#include <functional> 
#include <cctype>
#include <locale>
#include <string>
#include <vector>


#include "generalFunctions.h"


std::string& itostr(int n, std::string& s){
    if(n==0){
        s="0";
        return s;
    }
 
    int sign = -(n<0);
    unsigned int val = (n^sign)-sign;
 
    int size;
    if(val>=10000){
        if(val>=10000000){
            if(val>=1000000000)	size=10;
            else if(val>=100000000)	size=9;
            else size=8;
        }else{
            if(val>=1000000)	size=7;
            else if(val>=100000)	size=6;
            else	size=5;
        }
    }else{
        if(val>=100){
            if(val>=1000)	size=4;
            else	size=3;
        }else{
            if(val>=10) size=2;
            else size=1;
        }
    }
 
    s.resize(-sign+size);
    char* c = &s[0];
    if(sign) *c++='-';
 
    char* d = c+size-1;
    while(val>0){
        *d--='0' + (val % 10);
        val /= 10;
    }
    return s;
}

std::string int2str(int n){
	std::stringstream ss;
	ss << n;
	return ss.str();
}


std::string dbl2str(double n){
	std::stringstream ss;
	//ss.setf(std::ios::showpoint);
	ss << n;
	return ss.str();
}



std::string replace(const std::string old_string, const std::string new_string, std::string target_string) {
	for(std::string::size_type position = 0; (position = target_string.find(old_string, position)) != std::string::npos; ++position)
        target_string.replace(position, old_string.size(), new_string);
	return target_string;
}

// trim from start
std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

//// trim from both ends
std::string &bothtrim(std::string &s) {
        return ltrim(rtrim(s));
}

std::string trim(std::string str){
	std::string s = str;
	s.erase(s.find_last_not_of(" \n\r\t")+1);
	return s;
}


std::string toUpperCase(std::string str){
	transform(str.begin(), str.end(), str.begin(), toupper);
	return str;
}

std::vector<std::string> StringSplit(std::string strOrigin, std::string strTok, int limit){
	int cutAt = 0;
	int cnt = 0;
	std::vector<std::string>	strResult;

	if(limit == 0){
		while((cutAt = strOrigin.find_first_of(strTok)) != strOrigin.npos ){
			if(cutAt > 0)	strResult.push_back(strOrigin.substr(0,cutAt));
			strOrigin = strOrigin.substr(cutAt+1);
			cnt++;
		}
	}else{
		while((cutAt = strOrigin.find_first_of(strTok)) != strOrigin.npos || cnt < limit-1){
			if(cutAt > 0)	strResult.push_back(strOrigin.substr(0,cutAt));
			strOrigin = strOrigin.substr(cutAt+1);
			cnt++;
		}

	}

	if(strOrigin.length() > 0)	strResult.push_back(strOrigin.substr(0,cutAt));

	return strResult;
}






