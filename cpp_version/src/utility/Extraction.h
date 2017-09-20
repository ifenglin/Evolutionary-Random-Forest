//#define __USE_MINGW_ANSI_STDIO 0
//#define IMAGE_HEIGHT 1024
//#define IMAGE_WIDTH 2048
//#define IMAGE_CHANNELS 3
#include<complex>
#include<cstdlib>
#include<iostream>
#include<map>
#include<string>

int get_offset_regions(std::string, std::unique_ptr<int>&);
int make_features(std::string, std::unique_ptr<int>, int);