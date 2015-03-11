#include "qcsys.h"
#include <string>
extern void getPrefix(char*& pref);
const char* getIntcoFileName()
{
   static std::string strintco("");
   if (strintco.empty() ) {
      char* pref=NULL;
      getPrefix(pref);
      strintco = std::string(pref) + "intco.dat";
   }
   return strintco.c_str();
}

const char* getOptdataFileName()
{
   static std::string stroptdata("");
   if (stroptdata.empty() ) {
      char* pref=NULL;
      getPrefix(pref);
      stroptdata = std::string(pref) + "opt_data.1";
   }
   return stroptdata.c_str();
}

