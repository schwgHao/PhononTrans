#include <iostream>

extern "C" void vibrator_c(char* slabel, bool isbulk);

int main(){
	char* st="Cchain";
	vibrator_c(st, true);	
}
