#ifndef _PHYS_CONST_HH
#define _PHYS_CONST_HH

class PhysConst{
	public:
	PhysConst(std::string name,Double_t value,Double_t error,std::string units,std::string reference);  // constructor
	void SetReference(std::string reference);
	//~PhysConst(); // destructor (what for?)
	Double_t GetValue();   // get the value of the physical constant
	Double_t GetError();   // get the error of the physical constant
	std::string GetName(); 
	Double_t GetRandom();  // get a random value extracteg gauss around central value
	void PrintInfo(); // print to stderr a dump of the class object
	private:
	std::string name; // name of the constant
	Double_t value; // central value
	Double_t error; // uncertainty
	std::string units; // like [GeV]
	std::string reference;    // [OPTIONAL] the reference according to which the value was taken
};


#endif
