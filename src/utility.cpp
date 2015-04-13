#include "utility.h"
#define MAX_DATE 256

int get_date(void)
{
	SYSTEMTIME st;
	GetLocalTime(&st);
	std::string strMessage;

	char the_date[MAX_DATE];
	the_date[0] = '\0';

	sprintf(the_date,
		"%02d%02d%02d%02d%02d",
		//st.wYear,
		st.wMonth,
		st.wDay,
		st.wHour,
		st.wMinute,
		st.wSecond);

	strMessage = the_date;

	int date = (int)the_date;

	return  date;
}
double get_cpu_time(){
	FILETIME a, b, c, d;
	if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0){
		//  Returns total user time.
		//  Can be tweaked to include kernel times as well.
		return
			(double)(d.dwLowDateTime |
			((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
	}
	else{
		//  Handle error
		return 0;
	}
}



