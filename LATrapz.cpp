#include<iostream>
#include<vector>
#include<cmath>
#include<functional>

#include"LATrapz.hpp"
#include"misc.hpp"


using LD=double;
using vector=std::vector<LD>;

int main(){
	vector X=misc::linspace<LD>(0,1,100);
	std::function<LD(LD)> integrand=[](LD x){return x;};

	LATrapz<LD> trapz(integrand,X);
	
	std::cout << trapz.integrate() << "\n";
	
	
	
	
	return 0;
}