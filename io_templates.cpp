#pragma once 

inline void help_flag(string para,string help_string){
	// Padded string to 25 for format
	string full_help_string(40-para.length(), ' ');
	string para_str(para);
	full_help_string.insert(4,para_str);
	full_help_string.append(":");
	full_help_string.append(help_string);
	cout << full_help_string << endl;
} 
