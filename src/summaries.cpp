
#include "summaries.h"

void timesummary(time_t t1, time_t& t2)
{
	cout << "\nTotal time for computation:  ";

	if (t2>=3600)
	{
		cout << t2/3600 << " hours, ";
		t2 -= 3600*(t2/3600);
	}
	if (t2>=60)
	{
		cout << t2/60 << " minutes, ";
		t2 -= 60*(t2/60);
	} 
	if (t2<60)
		cout << t2 << (t2==1?" second":" seconds");
	if (t1>60)
		cout << " (" << t1 << " seconds).";
	cout << endl;
}

void errorMsg(char* filename)
{
	cerr << "Error: Could not open " << filename << "." << endl;
}
/*void errorMsg(char* suppText, char* filename) {
	cerr << "Error: could not open " << suppText << " " << filename << "." << endl;
}*/
