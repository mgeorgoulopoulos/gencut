#ifndef __EXCEPTION_H__
#define __EXCEPTION_H__

#include<string>
#include <vector>

class Exception {
public:
	using StringList = std::vector<std::string>;
	
	Exception() {}

	Exception(const std::string &message) : Exception(StringList{message}) {}
	
	Exception(const StringList &strings) : strings(strings) {}
	
	void print() {
		printf("Exception details:\n");
		for (const std::string &s : strings) {
			printf("\t%s\n", s.c_str());			
		}
	}

private:
	StringList strings;
	
};

#endif // __EXCEPTION_H__