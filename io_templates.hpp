#pragma once 

inline void help_flag(std::string para,std::string help_string){
	// Padded string to 25 for format
	std::string full_help_string(40-para.length(), ' ');
	std::string para_str(para);
	full_help_string.insert(4,para_str);
	full_help_string.append(":");
	full_help_string.append(help_string);
	std::cout << full_help_string << std::endl;
} 
