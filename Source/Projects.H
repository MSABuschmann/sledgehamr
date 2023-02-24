/* This file is auto-generated. MODIFICATIONS ARE FUTILE! */

#ifndef Projects_H_
#define Projects_H_

#include<AxionStrings.H>
#include<AxionStringsOther.H>

#define SLEDGEHAMR_PROJECT(str) {\
   if( str == "AxionStrings" ) return new AxionStrings;\
   if( str == "AxionStringsOther" ) return new AxionStringsOther;\
}

#endif
